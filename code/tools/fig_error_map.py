"""芯片误差地图(IBM error-map 风格):真实拓扑 + 双温度条(零额度)。

- 节点 = 比特,颜色 = **当日自测**读出误差 1-η(不是厂商标定——这是我们的方法论);
- 边 = 耦合器,颜色 = **自测 Bell 不保真** 1-bell_mean(全芯片图谱,仅 Baihua 有);
- 死比特/未自测比特 = 灰底红 ×;两条独立温度条(比特/边),各自单色渐变。

适用于全部三台机:Baihua(每日)、Dongling/Shenglian(有 e15 自测日;Shenglian
拓扑经网页同源端点 /sqc-api/chip-details 获取,见 qfab/client.py)。
η 来源:Baihua=e14 每日体检;其他机=results/e15_crossbackend/<机>_<日期>/ 自测。
拓扑快照缺失时自动在线拉取一次并落盘(零额度)。

产出 figures/<日期>/error_map_<机>.png/svg/pdf + 同名 nodes/edges CSV

运行(Mac): ./qrun tools/fig_error_map.py [Baihua|Dongling|...]
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize  # noqa: E402

from qfab import analysis  # noqa: E402
from tools import fig_style as st  # noqa: E402

CHIP = sys.argv[1] if len(sys.argv) > 1 else "Baihua"
VMAX_Q = 0.5    # 比特读出误差色标上限
VMAX_E = 0.4    # 边不保真色标上限

# 边用 teal 渐变(与节点蓝渐变区分开,两条温度条不打架)
_SEQ_TEAL = ["#E2F2F3", "#B7E0E3", "#8CCDD2", "#62BAC1", "#42949E",
             "#2F7A84", "#1F616A", "#124A52"]


def teal_cmap():
    return LinearSegmentedColormap.from_list("seq_teal", _SEQ_TEAL)


def latest_snapshot() -> dict:
    files = sorted((PROJ / "results" / "chip_snapshots").glob(f"chip_{CHIP}_*.json"))
    if files:
        return json.loads(files[-1].read_text())
    from qfab import client                # 无缓存 → 在线拉一次(零额度,自动落盘)
    return client.get_chip_info(CHIP, refresh=True)


def latest_bell() -> dict[frozenset, float]:
    files = sorted((PROJ / "results" / "e03_bell").glob(f"atlas_{CHIP}_*.json"))
    if not files:
        return {}
    rows = json.loads(files[-1].read_text()).get("rows", [])
    return {frozenset(r["edge"]): r["bell_mean"] for r in rows
            if r.get("bell_mean") is not None}


def latest_eta() -> dict:
    """{date, time, eta:{q:η}}:Baihua=e14 体检;其他机=e15 跨机自测。
    e15 自校准可能按宽度切成多组(selfcal_prep0_g0/_g1/...),按 meta.qubits
    配对 prep0/prep1 并合并所有组的 η。"""
    if CHIP == "Baihua":
        from tools.fig_weather_chip import frames
        return frames()[-1]
    for d in sorted((PROJ / "results" / "e15_crossbackend").glob(f"{CHIP}_*"),
                    reverse=True):
        pairs: dict[tuple, dict[int, dict]] = {}
        for f in sorted(d.glob("selfcal_prep*__*.json")):
            r = json.loads(f.read_text())
            if not r.get("result", {}).get("count"):
                continue
            key = tuple(r["meta"]["qubits"])
            prep = 1 if r["meta"]["prep"] == "prep1" else 0
            pairs.setdefault(key, {})[prep] = r
        eta: dict[int, float] = {}
        for qubits, pr in pairs.items():
            if 0 not in pr or 1 not in pr:
                continue
            asg = analysis.assignment_from_p1(
                analysis.bit_p1(pr[0]["result"]["count"], len(qubits)),
                analysis.bit_p1(pr[1]["result"]["count"], len(qubits)))
            eta.update({q: asg[i]["eta"] for i, q in enumerate(qubits)})
        if eta:
            return {"date": d.name.split("_", 1)[1], "time": "", "eta": eta}
    sys.exit(f"{CHIP} 无自测 η 数据(先跑 e15 跨机自校准)")


def main():
    st.apply_style()
    snap = latest_snapshot()
    fr = latest_eta()
    bell = latest_bell()

    cols = snap["size"][1]
    pos = {q: (q % cols, -(q // cols)) for q in range(len(snap["qubits_info"]))}

    fig, ax = plt.subplots(figsize=(9.2, 7.4))
    ax.grid(False)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    # ---- 边
    n_q = Normalize(0, VMAX_Q)
    n_e = Normalize(0, VMAX_E)
    ecmap = teal_cmap()
    # 口径纪律:边的颜色只用自测 Bell 不保真;厂商 CZ 是另一种口径(门保真,
    # 数值低一个量级),混同一条温度条会误导 → 只进 CSV,不上色
    segs_meas, cols_meas, segs_un = [], [], []
    edge_rows = []
    for c in snap["couplers_info"].values():
        a, b = c["qubits_index"]
        seg = [pos[a], pos[b]]
        key = frozenset((a, b))
        if key in bell:
            infid = 1 - bell[key]
            segs_meas.append(seg); cols_meas.append(ecmap(n_e(min(infid, VMAX_E))))
        else:
            infid = None
            segs_un.append(seg)
        edge_rows.append({"q_a": a, "q_b": b,
                          "bell_infidelity_selfcal": infid,
                          "vendor_cz_fidelity": c.get("fidelity")})
    ax.add_collection(LineCollection(segs_un, colors=st.ROLE["base2"],
                                     linewidths=1.6, zorder=1,
                                     linestyle=(0, (2, 1.5))))
    ax.add_collection(LineCollection(segs_meas, colors=cols_meas, linewidths=4.2,
                                     zorder=2))

    # ---- 节点
    qcmap = st.seq_cmap()
    node_rows = []
    for q, (x, y) in pos.items():
        if q in fr["eta"]:
            err = 1 - fr["eta"][q]
            ax.scatter([x], [y], s=190, color=qcmap(n_q(min(err, VMAX_Q))),
                       zorder=3, edgecolors=st.SURFACE, linewidths=1.4)
            ink = "#ffffff" if err > 0.5 * VMAX_Q else st.INK2
            ax.text(x, y, str(q), ha="center", va="center", fontsize=5.2,
                    color=ink, zorder=4)
            node_rows.append({"qubit": q, "readout_err": round(err, 4),
                              "eta": round(fr["eta"][q], 4), "status": "ok"})
        else:
            ax.scatter([x], [y], s=190, color=st.ROLE["base2"], zorder=3,
                       edgecolors=st.SURFACE, linewidths=1.4)
            ax.text(x, y, "×", ha="center", va="center", fontsize=8,
                    color=st.CRITICAL, fontweight="bold", zorder=4)
            node_rows.append({"qubit": q, "readout_err": None, "eta": None,
                              "status": "unavailable"})
    ax.set_aspect("equal")
    ax.set_xlim(-0.7, cols - 0.3)
    ax.set_ylim(-(len(pos) // cols) + 0.3, 0.7)

    cb1 = fig.colorbar(ScalarMappable(norm=n_q, cmap=qcmap), ax=ax,
                       shrink=0.72, pad=0.02)
    cb1.set_label("比特读出误差 1-η(当日自测)", color=st.INK2)
    cb1.outline.set_visible(False)
    cb2 = fig.colorbar(ScalarMappable(norm=n_e, cmap=ecmap), ax=ax,
                       shrink=0.72, pad=0.06)
    cb2.set_label("双比特边 Bell 不保真(当日自测;灰虚线=当日未自测)",
                  color=st.INK2)
    cb2.outline.set_visible(False)

    # 灰边成分注脚:未自测 ≠ 没测,而是平台不可用(触死比特/厂商未标定→拒收)
    dead_v = {v["index"] for v in snap["qubits_info"].values()
              if not v.get("T1") or not v.get("fidelity")}
    n_deadedge = sum(1 for c in snap["couplers_info"].values()
                     if set(c["qubits_index"]) & dead_v)
    n_uncal = sum(1 for c in snap["couplers_info"].values()
                  if not (set(c["qubits_index"]) & dead_v) and not c.get("fidelity"))
    n_all = len(snap["couplers_info"])
    if segs_meas:
        cap = (f"边自测 {len(segs_meas)}/{n_all};灰虚线 {len(segs_un)} 条中 "
               f"{n_deadedge} 条触不可用比特、{n_uncal} 条厂商未标定(平台拒收,"
               f"即当日全部可用边均已自测)")
    else:
        cap = f"当日无该机边自测数据(灰虚线={n_all} 条边;跑 e03 边图谱后自动上色)"
    ax.text(0, -0.05, cap, transform=ax.transAxes, fontsize=8.5, color=st.MUTED)

    ax.set_title(f"芯片误差地图 · {CHIP} · 自测 {fr['date']} {fr['time']} · "
                 f"标定时刻 {snap.get('calibration_time', '?')}", pad=12)
    st.save(fig, f"error_map_{CHIP}", formats=("png", "svg", "pdf"))

    for name, rows_, cols_ in ((f"error_map_{CHIP}_nodes", node_rows,
                                ["qubit", "readout_err", "eta", "status"]),
                               (f"error_map_{CHIP}_edges", edge_rows,
                                ["q_a", "q_b", "bell_infidelity_selfcal",
                                 "vendor_cz_fidelity"])):
        p = st.outdir() / f"{name}.csv"
        with p.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols_)
            w.writeheader()
            w.writerows(rows_)
        print(f"已写 {p}")


if __name__ == "__main__":
    main()
