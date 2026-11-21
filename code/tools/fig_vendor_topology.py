"""官方拓扑图(零额度):厂商标定视图,与我们的自测误差地图并排
= "快照 vs 实测"方法论的图版。Baihua 上同时高亮我们当日选址的 d=7 链(官方图上
看我们的锚点落位)。

Baihua/Dongling 用平台 SDK 自带渲染(quark Backend.draw);Shenglian 在 SDK
白名单之外(quarkcircuit 0.5.5 只认 Baihua/Dongling/Miaofeng)→ 用 chip-details
端点的官方标定 JSON(含比特网格坐标)自渲染,数据口径同为厂商标定。

产出 figures/<日期>/vendor_topology_<机>.svg(+ Shenglian 版配同名 CSV)

运行(联网,Mac/服务器均可): ./qrun tools/fig_vendor_topology.py
"""
from __future__ import annotations

import sys
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))

from qfab import client  # noqa: E402
from tools import fig_style as st  # noqa: E402


def main():
    client.get_task()
    # 注:不走 quark.Task.backend()——其包装层按位置传参,highlight_nodes/
    # save_svg_fname 会错位塞进 show_*_index(平台 SDK bug,已记反馈素材);
    # 直接关键字调用 Backend.draw 绕过。
    from quark.circuit import Backend
    st.apply_style()
    out = st.outdir()

    jobs = []
    for chip in ("Baihua", "Dongling"):
        hl = []
        if chip == "Baihua":
            try:
                from qfab import anchors as anchmod
                hl = list(anchmod.latest(chip).get("chains", {}).get("7", []))
            except Exception:  # noqa: BLE001
                hl = []
        jobs.append((chip, hl))

    import matplotlib.pyplot as plt

    def render(chip, hl, suffix):
        f = out / f"vendor_topology_{chip}{suffix}"
        try:
            Backend(chip).draw(show_couplers_fidelity=True,
                               show_qubits_attributes="fidelity",
                               show_qubits_index=True,
                               highlight_nodes=hl,
                               save_svg_fname=str(f))
            plt.close("all")
            svg = f.with_suffix(".svg")
            if svg.exists():
                print(f"已写 {svg}({svg.stat().st_size // 1024} KB)"
                      + (f",高亮当日链 {hl}" if hl else "(全图)"))
                import subprocess
                subprocess.run(["qlmanage", "-t", "-s", "2400", "-o",
                                str(out), str(svg)],
                               capture_output=True)  # macOS 缩略图 → 同名 .svg.png
            else:
                print(f"⚠️ {chip}: 调用成功但未见 {svg}")
        except Exception as e:  # noqa: BLE001
            print(f"⚠️ {chip} 官方拓扑图失败: {str(e)[:100]}")

    for chip, hl in jobs:
        render(chip, [], "")                 # 全图版(=官方 demo 同款)
        if hl:
            render(chip, hl, "_chain")       # 高亮版:当日选址链落在官方地图上
    render_shenglian(out)


def render_shenglian(out):
    """Shenglian:SDK 白名单外 → 用 chip-details 官方标定自渲染(厂商口径)。"""
    import csv
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    from qfab import client as cl

    info = cl.get_chip_info("Shenglian", refresh=True)
    q = {v["index"]: v for v in info["qubits_info"].values()}
    pos = {i: (v["coordinate"][0], -v["coordinate"][1]) for i, v in q.items()}

    fig, ax = plt.subplots(figsize=(9.6, 8.2))
    ax.grid(False); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    n_e = Normalize(0.95, 1.0)
    ecmap = st.seq_cmap()
    edge_rows = []
    for c in info["couplers_info"].values():
        a, b = c["qubits_index"]
        f = c.get("fidelity") or 0
        (x1, y1), (x2, y2) = pos[a], pos[b]
        if f > 0:
            ax.plot([x1, x2], [y1, y2], color=ecmap(n_e(max(f, 0.95))),
                    lw=4.0, zorder=2, solid_capstyle="round")
        else:
            ax.plot([x1, x2], [y1, y2], color=st.ROLE["base2"], lw=1.4,
                    ls=(0, (2, 1.5)), zorder=1)
        edge_rows.append({"q_a": a, "q_b": b, "vendor_cz_fidelity": f or None})

    node_rows = []
    for i, v in q.items():
        x, y = pos[i]
        ro = ((v.get("readout g_fidelity") or 0)
              + (v.get("readout e_fidelity") or 0)) / 2
        dead = (v.get("T1") or 0) <= 0 or (v.get("fidelity") or 0) <= 0
        face = st.ROLE["base2"] if (dead or ro == 0) else st.SURFACE
        ax.scatter([x], [y], s=210, color=face, zorder=3,
                   edgecolors=st.INK2, linewidths=0.9)
        ax.text(x, y, "×" if dead else str(i), ha="center", va="center",
                fontsize=8 if dead else 5.4,
                color=st.CRITICAL if dead else st.INK2,
                fontweight="bold" if dead else "normal", zorder=4)
        node_rows.append({"qubit": i, "T1_us": v.get("T1"), "T2_us": v.get("T2"),
                          "gate_fidelity_1q": v.get("fidelity"),
                          "readout_mean_fidelity": round(ro, 4) or None,
                          "status": "dead" if dead else
                                    ("readout_uncal" if ro == 0 else "ok")})
    ax.set_aspect("equal")
    cb = fig.colorbar(ScalarMappable(norm=n_e, cmap=ecmap), ax=ax,
                      shrink=0.7, pad=0.02)
    cb.set_label("厂商 CZ 门保真度(官方标定)", color=st.INK2)
    cb.outline.set_visible(False)
    n_cal = sum(1 for r in edge_rows if r["vendor_cz_fidelity"])
    ax.text(0, -0.04,
            f"已标定 {n_cal}/{len(edge_rows)} 条边;灰虚线=厂商未标定;"
            "灰底节点=读出未标定;×=不可用。数据=官方 chip-details 标定"
            "(与 Dongling 同构的 12×7 对角格)",
            transform=ax.transAxes, fontsize=8.5, color=st.MUTED)
    ax.set_title(f"官方拓扑与标定 · Shenglian · "
                 f"标定时刻 {info.get('calibration_time', '?')}", pad=12)
    st.save(fig, "vendor_topology_Shenglian", formats=("png", "svg", "pdf"))
    for name, rows_, cols_ in (("vendor_topology_Shenglian_nodes", node_rows,
                                ["qubit", "T1_us", "T2_us", "gate_fidelity_1q",
                                 "readout_mean_fidelity", "status"]),
                               ("vendor_topology_Shenglian_edges", edge_rows,
                                ["q_a", "q_b", "vendor_cz_fidelity"])):
        p = out / f"{name}.csv"
        with p.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols_)
            w.writeheader(); w.writerows(rows_)
        print(f"已写 {p}")


if __name__ == "__main__":
    main()
