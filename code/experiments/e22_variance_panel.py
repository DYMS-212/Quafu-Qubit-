"""e22 · 编码 = 方差抑制器:固定面板上量 σ(编码) vs σ(裸)。

发现(e19):编码逻辑比特错误率稳(~0.05),裸物理比特 volatile(0.02-0.09),门级
"交叉"与否主要取决于裸比特那天多差。要把这升级成硬结论"编码买的是可预测性",必须
**固定一组跨日不变的 star+bare 面板**每天跑(否则比特身份漂移与器件漂移纠缠)。

本实验:
  - 固定面板(首次生成后冻结到 results/anchors/e22_panel.json,永不变)
  - 每轮跑 [[4,2,2]] 编码 + 同面板裸物理比特,各 k∈{0,4,8} 层 transversal H
  - 多轮/多天 pool_ci → σ(enc) vs σ(bare):σ(enc) ≪ σ(bare) 即"编码=可预测性"
运行: ./qrun experiments/e22_variance_panel.py
"""
from __future__ import annotations

import glob
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, chip as chipmod, client, provenance  # noqa: E402
from qfab.qasm import QASM  # noqa: E402

CHIP = "Baihua"
RUN = time.strftime("%Y-%m-%d_%H%M")
EXP = f"e22_variance_panel/{RUN}"
PROJ = Path(__file__).resolve().parent.parent
OUT = PROJ / "results" / "e22_variance_panel"
PANEL_FILE = PROJ / "results" / "anchors" / "e22_panel.json"
SHOTS = 4096
KS = [0, 4, 8]


def load_today_eta() -> dict[int, float]:
    d = PROJ / "results" / "e14_probe" / time.strftime("%Y-%m-%d")
    rec0 = rec1 = None
    for f in d.glob("selfcal_prep0__*.json"):
        rec0 = json.loads(f.read_text())
    for f in d.glob("selfcal_prep1__*.json"):
        rec1 = json.loads(f.read_text())
    if not (rec0 and rec1):
        raise RuntimeError("缺当日自校准,先跑 e14_daily_probe")
    qubits = rec0["meta"]["qubits"]
    asg = analysis.assignment_from_p1(
        analysis.bit_p1(rec0["result"]["count"], len(qubits)),
        analysis.bit_p1(rec1["result"]["count"], len(qubits)))
    return {q: asg[i]["eta"] for i, q in enumerate(qubits)}


def get_panel(ch, eta) -> dict:
    """固定面板:首次用当日 η 选一个好星 + 星外一个好裸比特,冻结;之后永远读它。"""
    if PANEL_FILE.exists():
        return json.loads(PANEL_FILE.read_text())
    av = set(ch.available_qubits()) & set(eta)
    best = None
    for c in av:
        nb = [x for x in ch.neighbors(c) if x in av and (ch.edge_fid(c, x) or 0) > 0.5]
        if len(nb) >= 3:
            leaves = sorted(nb, key=lambda x: -eta[x])[:3]
            s = min(eta[q] for q in (c, *leaves))
            if best is None or s > best[0]:
                best = (s, c, leaves)
    star = [best[1], *best[2]]
    bare = next(q for q in sorted(av, key=lambda q: -eta[q]) if q not in star)
    panel = {"star": star, "bare": [bare], "frozen_at": time.strftime("%Y-%m-%d %H:%M"),
             "chip_calibration_time": ch.calibration_time,
             "note": "e22 固定方差面板,永不修改(与冻结锚点同理)"}
    PANEL_FILE.parent.mkdir(parents=True, exist_ok=True)
    PANEL_FILE.write_text(json.dumps(panel, ensure_ascii=False, indent=2))
    print(f"⚠️ e22 面板首次生成并冻结 → {PANEL_FILE}")
    return panel


def circ_encoded(k: int) -> str:
    qc = QASM(4)
    qc.h(0).cx(0, 1).cx(0, 2).cx(0, 3)
    qc.barrier()
    for _ in range(k):
        for q in range(4):
            qc.h(q)
        qc.barrier()
    if k % 2 == 1:
        for q in range(4):
            qc.h(q)
    qc.barrier().measure_all()
    return qc.qasm()


def circ_bare(k: int) -> str:
    qc = QASM(1)
    qc.barrier()
    for _ in range(k):
        qc.h(0)
        qc.barrier()
    if k % 2 == 1:
        qc.h(0)
    qc.barrier().measure_all()
    return qc.qasm()


def eval_encoded(count: dict, k: int) -> dict:
    tot = acc = e = 0
    for key, v in count.items():
        b = analysis.bits_of(key, 4)
        tot += v
        if (b[0] ^ b[1] ^ b[2] ^ b[3]) == 0:
            acc += v
            l1, l2 = (b[0] ^ b[2], b[0] ^ b[1]) if k % 2 == 0 else (b[0] ^ b[1], b[0] ^ b[2])
            e += v * ((l1 != 0) + (l2 != 0))
    return {"acceptance": round(acc / (tot or 1), 4),
            "err": round(e / (2 * (acc or 1)), 5),
            "err_ci": [round(x, 5) for x in analysis.wilson_ci(e, 2 * (acc or 1))]}


def main():
    print(f"后端状态: {client.status()}")
    ch = chipmod.load(CHIP, refresh=True)
    panel = get_panel(ch, load_today_eta())
    star, bare = panel["star"], panel["bare"]
    print(f"固定面板: star={star} bare={bare}(冻结于 {panel.get('frozen_at')})")

    specs = []
    for k in KS:
        specs.append(dict(circuit=circ_encoded(k), chip=CHIP, shots=SHOTS,
                          name=f"e22_enc_k{k}", exp=EXP, label=f"enc_k{k}",
                          options={"target_qubits": star},
                          meta={"arm": "enc", "k": k}))
        specs.append(dict(circuit=circ_bare(k), chip=CHIP, shots=SHOTS,
                          name=f"e22_bare_k{k}", exp=EXP, label=f"bare_k{k}",
                          options={"target_qubits": bare},
                          meta={"arm": "bare", "k": k}))
    recs = client.run_batch(specs, timeout=1200)

    rows = []
    for r in recs:
        c = client.counts(r)
        if not c:
            continue
        m = r["meta"]
        if m["arm"] == "enc":
            ev = eval_encoded(c, m["k"])
        else:
            p1 = analysis.bit_p1(c, 1)[0]
            sh = analysis.shots_of(c)
            ev = {"err": round(p1, 5),
                  "err_ci": [round(x, 5) for x in analysis.wilson_ci(round(p1 * sh), sh)]}
        rows.append({"arm": m["arm"], "k": m["k"], "tid": r.get("tid"), **ev})

    OUT.mkdir(parents=True, exist_ok=True)
    summary = {"provenance": provenance.stamp(
                   anchor=None, chip=CHIP, tids=provenance.tids_of(recs),
                   qubits={"star": star, "bare": bare},
                   extra={"panel_frozen_at": panel.get("frozen_at"),
                          "chip_calibration_time": ch.calibration_time}),
               "rows": rows}
    (OUT / f"summary_{RUN}.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2))

    # 跨轮/跨日聚合:σ(enc) vs σ(bare)(读全部 e22 summary)
    print("\n===== 固定面板:σ(编码) vs σ(裸)——跨全部轮次 =====")
    hist = {("enc", k): [] for k in KS}
    hist.update({("bare", k): [] for k in KS})
    for f in glob.glob(str(OUT / "*" / "summary_*.json")) + glob.glob(str(OUT / "summary_*.json")):
        try:
            d = json.loads(Path(f).read_text())
        except Exception:  # noqa: BLE001
            continue
        for rr in d.get("rows", []):
            key = (rr.get("arm"), rr.get("k"))
            if key in hist and rr.get("err") is not None:
                hist[key].append(rr["err"])
    for k in KS:
        pe, pb = analysis.pool_ci(hist[("enc", k)]), analysis.pool_ci(hist[("bare", k)])
        if pe.get("n") and pb.get("n"):
            print(f"  k={k}: enc σ={pe.get('sd', 0):.4f}(n={pe['n']}, {pe['best']:.4f}-{pe['worst']:.4f}) | "
                  f"bare σ={pb.get('sd', 0):.4f}(n={pb['n']}, {pb['best']:.4f}-{pb['worst']:.4f})")
    print(f"已写 {OUT / f'summary_{RUN}.json'}")


if __name__ == "__main__":
    main()
