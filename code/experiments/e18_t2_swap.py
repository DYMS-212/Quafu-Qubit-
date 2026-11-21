"""e18 · T=2 多轮 syndrome 的 SWAP 合成(时间盒项:全场唯一"跨时间纠错"的机会)。

背景:平台无 mid-circuit measurement;且图论上(heavy-hex 公共邻居 max=1)
native fresh-ancilla T=2 不存在。唯一窄路:round1 的 ancilla 用 SWAP(3 CX)
换出到旁边的 flag 位保存,原位 ancilla 复用做 round2,末端一起测量。

线路(q0=D0, q1=D1, q2=A, q3=F,F 是 A 的第三邻居):
  round1: CX(D0→A) CX(D1→A)
  swap:   SWAP(A,F)                  ← round1 syndrome 存进 F
  [两轮之间的时间窗:mid 注错在这里]
  round2: CX(D0→A) CX(D1→A)
  测量:  F=round1 syndrome,A=round2 syndrome

查找表(期望):
  none      → (s1,s2)=(0,0)
  X_D0_pre  → (1,1)   编码前的错,两轮都看见
  X_D0_mid  → (0,1)   两轮之间的错,只有 round2 看见 ← 时间定位!
  X_A_mid   → (0,0)+A读出翻转 → 对照
成功标准:mid 错的 (0,1) 特征显著高于随机 → 跨时间信息真实存在。

运行: ./qrun experiments/e18_t2_swap.py
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
RUN = time.strftime("%Y-%m-%d_%H%M")
EXP = f"e18_t2swap/{RUN}"
PROJ = Path(__file__).resolve().parent.parent
OUT = PROJ / "results" / "e18_t2swap"
SHOTS = 8192


def find_layout(ch, eta):
    """找 D0-A-D1 三元组 + A 的第三邻居 F,按 min-η 打分。"""
    av = set(ch.available_qubits()) & set(eta)
    best = None
    for a in av:
        nb = [x for x in ch.neighbors(a) if x in av
              and (ch.edge_fid(a, x) or 0) > 0.5]
        if len(nb) >= 3:
            nb = sorted(nb, key=lambda x: -eta[x])
            d0, d1, f = nb[0], nb[1], nb[2]
            score = min(eta[q] for q in (a, d0, d1, f))
            if best is None or score > best[0]:
                best = (score, d0, d1, a, f)
    return best[1:] if best else None


def circ(case: str) -> str:
    qc = QASM(4)  # q0=D0 q1=D1 q2=A q3=F
    if case == "X_D0_pre":
        qc.x(0)
    qc.barrier()
    qc.cx(0, 2).cx(1, 2)          # round 1
    qc.barrier()
    qc.cx(2, 3).cx(3, 2).cx(2, 3)  # SWAP(A,F)
    qc.barrier()
    if case == "X_D0_mid":
        qc.x(0)
    if case == "X_D1_mid":
        qc.x(1)
    if case == "X_A_mid":
        qc.x(2)
    qc.barrier()
    qc.cx(0, 2).cx(1, 2)          # round 2
    qc.barrier().measure_all()
    return qc.qasm()


def main():
    print(f"后端状态: {client.status()}")
    files = sorted((PROJ / "results" / "e02_selfcal").glob("assignment_*.json"))
    eta = {r["qubit"]: r["eta"] for r in json.loads(files[-1].read_text())}
    ch = chipmod.load(CHIP, refresh=True)
    layout = find_layout(ch, eta)
    if not layout:
        print("找不到 D0-A-D1+F 布局")
        return
    d0, d1, a, f = layout
    tq = [d0, d1, a, f]
    print(f"布局: D0=Q{d0} D1=Q{d1} A=Q{a} F=Q{f} "
          f"(η={[round(eta[q], 3) for q in tq]})")

    cases = ("none", "X_D0_pre", "X_D0_mid", "X_D1_mid", "X_A_mid")
    specs = [dict(circuit=circ(cs), chip=CHIP, shots=SHOTS,
                  name=f"t2_{cs}", exp=EXP, label=cs,
                  options={"target_qubits": tq},
                  meta={"case": cs, "layout": tq}) for cs in cases]
    recs = client.run_batch(specs, timeout=1200)

    # X_A_mid:A 的 round1 值已 SWAP 到 F,翻转只污染 round2 → 期望 (0,1)。
    # 注意与 X_D0_mid 特征简并(单轮内不可分辨,完整多轮才能消歧)——诚实局限,写进报告。
    expect = {"none": (0, 0), "X_D0_pre": (1, 1),
              "X_D0_mid": (0, 1), "X_D1_mid": (0, 1), "X_A_mid": (0, 1)}
    rows = []
    print("\n===== T=2 时空查找表(s1=F, s2=A)=====")
    for r in recs:
        c = client.counts(r)
        if not c:
            continue
        cs = r["meta"]["case"]
        tot = analysis.shots_of(c)
        joint = {}
        hit = 0
        for key, v in c.items():
            b = analysis.bits_of(key, 4)
            s1, s2 = b[3], b[2]
            joint[f"{s1}{s2}"] = joint.get(f"{s1}{s2}", 0) + v
            if (s1, s2) == expect[cs]:
                hit += v
        dist = {k: round(v / tot, 4) for k, v in sorted(joint.items())}
        rows.append({"case": cs, "expected": list(expect[cs]),
                     "p_expected": round(hit / tot, 4),
                     "p_expected_ci": [round(x, 4) for x in analysis.wilson_ci(hit, tot or 1)],
                     "tid": r.get("tid"), "dist": dist})
        print(f"  {cs:9s} 期望{expect[cs]}: P={hit / tot:.4f}  分布={dist}")
    OUT.mkdir(parents=True, exist_ok=True)
    summary = {"provenance": provenance.stamp(
                   anchor=None, qubits=tq, tids=provenance.tids_of(recs), chip=CHIP,
                   extra={"anchor": "e02_selfcal_eta + find_layout",
                          "chip_calibration_time": ch.calibration_time}),
               "layout": tq, "rows": rows}
    (OUT / f"summary_{RUN}.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"已写 {OUT / f'summary_{RUN}.json'}")


if __name__ == "__main__":
    main()
