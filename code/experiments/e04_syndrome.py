"""e04 · syndrome 原语:weight-3 星形 Z/X 检验注错查找表(heavy-hex 原生最大检验)。

Z 检验:叶子=data(|0>),中心=ancilla:CX(leaf→anc)×3;注入 X 于单个 data。
X 检验:叶子=data(|+>),中心 H·CX(anc→leaf)×3·H(相位反冲测 XXX);注入 Z。
(weight-2 ZZ 三元组已由 e14 每日体检覆盖,此处不重复。)

星形选址:e02 最新自测 η 打分,取最好且边 vendor CZ>0.5 的一颗。

运行: ./qrun experiments/e04_syndrome.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, chip as chipmod, client, provenance  # noqa: E402
from qfab.qasm import QASM  # noqa: E402

CHIP = "Baihua"
TODAY = time.strftime("%Y-%m-%d")
RUN = time.strftime("%Y-%m-%d_%H%M")   # 同日多轮互不覆盖
EXP = f"e04_syndrome/{RUN}"
PROJ = Path(__file__).resolve().parent.parent
OUT = PROJ / "results" / "e04_syndrome"
SHOTS = 4096


def load_eta() -> tuple[dict[int, float], str]:
    files = sorted((PROJ / "results" / "e02_selfcal").glob("assignment_*.json"))
    rows = json.loads(files[-1].read_text())
    return {r["qubit"]: r["eta"] for r in rows}, files[-1].name


def best_star(ch, eta):
    av = set(ch.available_qubits()) & set(eta)
    best = None
    for c in av:
        nb = [x for x in ch.neighbors(c) if x in av
              and (ch.edge_fid(c, x) or 0) > 0.5]
        if len(nb) >= 3:
            leaves = sorted(nb, key=lambda x: -eta[x])[:3]
            score = min(eta[q] for q in [c, *leaves])
            if best is None or score > best[0]:
                best = (score, c, leaves)
    return best[1], best[2]


def circ_z(inject: int | None) -> str:
    """q0..q2=data(叶),q3=ancilla(中心)。Z 检验:anc = ⊕data。"""
    qc = QASM(4)
    if inject is not None:
        qc.x(inject)
    qc.barrier()
    for i in range(3):
        qc.cx(i, 3)
    qc.barrier().measure_all()
    return qc.qasm()


def circ_x(inject: int | None) -> str:
    """X 检验:data 制备 |+>,相位反冲到 ancilla;注入 Z 翻转对应 X 本征值。"""
    qc = QASM(4)
    for i in range(3):
        qc.h(i)
    qc.barrier()
    if inject is not None:
        qc.z(inject)
    qc.barrier()
    qc.h(3)
    for i in range(3):
        qc.cx(3, i)
    qc.h(3)
    qc.barrier().measure_all()
    return qc.qasm()


def main():
    print(f"后端状态: {client.status()}")
    ch = chipmod.load(CHIP, refresh=True)
    eta, eta_src = load_eta()
    center, leaves = best_star(ch, eta)
    tq = [*leaves, center]
    print(f"星形: 中心 Q{center} + 叶 {leaves} "
          f"(η={[round(eta[q], 3) for q in tq]})")

    specs = []
    for check, builder in (("Z", circ_z), ("X", circ_x)):
        for inj in (None, 0, 1, 2):
            tag = "none" if inj is None else f"inj{inj}"
            specs.append(dict(circuit=builder(inj), chip=CHIP, shots=SHOTS,
                              name=f"syn3_{check}_{tag}", exp=EXP,
                              label=f"{check}_{tag}",
                              options={"target_qubits": tq},
                              meta={"check": check, "inject": inj,
                                    "center": center, "leaves": leaves}))
    recs = client.run_batch(specs, timeout=1200)

    rows = []
    print("\n===== weight-3 星形检验查找表(anc 期望: none→0, 单错→1)=====")
    for r in recs:
        c = client.counts(r)
        if not c:
            continue
        m = r["meta"]
        p1_anc = analysis.bit_p1(c, 4)[3]
        expected = 0 if m["inject"] is None else 1
        success = (1 - p1_anc) if expected == 0 else p1_anc
        tot = analysis.shots_of(c)
        s_lo, s_hi = analysis.wilson_ci(round(success * tot), tot)
        a_lo, a_hi = analysis.wilson_ci(round(p1_anc * tot), tot)
        rows.append({**{k: m[k] for k in ("check", "inject")},
                     "p_anc_1": round(p1_anc, 4),
                     "p_anc_1_ci": [round(a_lo, 4), round(a_hi, 4)],
                     "success": round(success, 4),
                     "success_ci": [round(s_lo, 4), round(s_hi, 4)],
                     "shots": tot})
        tag = "none" if m["inject"] is None else f"inj{m['inject']}"
        print(f"  {m['check']} {tag:5s}: P(anc=期望)={success:.4f}")
    prov = provenance.stamp(
        qubits=tq, tids=provenance.tids_of(recs), chip=CHIP,
        extra={"selection": "e02_eta_best_star (min-η, edge CZ>0.5)",
               "eta_source": eta_src, "center": center, "leaves": leaves,
               "shots_per_task": SHOTS,
               "ci": "Wilson 95% on P(anc=expected) 与 P(anc=1)"})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"summary_{RUN}.json").write_text(
        json.dumps({"provenance": prov, "center": center, "leaves": leaves,
                    "rows": rows}, ensure_ascii=False, indent=2))
    print(f"已写 {OUT / f'summary_{RUN}.json'}")


if __name__ == "__main__":
    main()
