"""e13c · 读出缓解的真机验证(held-out 图案恢复测试)。

缓解矩阵由全|0>/全|1>标定(e14)——但真实读出可能依赖上下文(邻居状态)。
本实验制备标定集之外的已知图案(随机比特串×4 + 棋盘格×2)于当日最优 d7 链
(13 比特),对比 raw / 缓解后 的图案恢复保真度:
  恢复到位 → 张量缓解模型在真机成立;残差 → 上下文效应的实测偏差(诚实边界)。

运行: ./qrun experiments/e13c_mitigation_validation.py
"""
from __future__ import annotations

import json
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, anchors as anchmod, client, mitigation, provenance  # noqa: E402
from qfab.qasm import QASM  # noqa: E402

CHIP = "Baihua"
RUN = time.strftime("%Y-%m-%d_%H%M")
EXP = f"e13_mitigation/val_{RUN}"
PROJ = Path(__file__).resolve().parent.parent
OUT = PROJ / "results" / "e13_mitigation"
SHOTS = 4096


def load_cal() -> dict[int, dict]:
    files = sorted((PROJ / "results" / "e02_selfcal").glob("assignment_*.json"))
    return {r["qubit"]: r for r in json.loads(files[-1].read_text())}


def pattern_circ(bits: list[int]) -> str:
    qc = QASM(len(bits))
    for i, b in enumerate(bits):
        if b:
            qc.x(i)
    qc.barrier().measure_all()
    return qc.qasm()


def main():
    print(f"后端状态: {client.status()}")
    a = anchmod.latest(CHIP)
    reg = a["chains"]["7"]              # 当日最优 13 比特链做验证寄存器
    n = len(reg)
    cal = load_cal()
    cal_src = sorted((PROJ / "results" / "e02_selfcal").glob("assignment_*.json"))[-1].name
    rng = random.Random(20260703)
    patterns = {f"rand{i}": [rng.randint(0, 1) for _ in range(n)] for i in range(4)}
    patterns["checker"] = [i % 2 for i in range(n)]
    patterns["checker_inv"] = [(i + 1) % 2 for i in range(n)]

    specs = [dict(circuit=pattern_circ(bits), chip=CHIP, shots=SHOTS,
                  name=f"mval_{tag}", exp=EXP, label=tag,
                  options={"target_qubits": reg},
                  meta={"pattern": bits, "qubits": reg})
             for tag, bits in patterns.items()]
    recs = client.run_batch(specs, timeout=1200)

    fs = mitigation.fs_for_qubits(reg, cal)
    rows = []
    print(f"\n===== held-out 图案恢复(寄存器={reg})=====")
    for r in recs:
        c = client.counts(r)
        if not c:
            continue
        bits = r["meta"]["pattern"]
        target = "".join(str(b) for b in reversed(bits))  # c_{n-1}..c0
        cm, clip = mitigation.correct_counts(c, fs)

        tot = analysis.shots_of(c)
        k_raw = sum(v for kk, v in c.items() if kk.zfill(n) == target)
        f_raw = k_raw / (tot or 1)
        f_mit = (sum(v for kk, v in cm.items() if kk.zfill(n) == target)
                 / (analysis.shots_of(cm) or 1))
        r_lo, r_hi = analysis.wilson_ci(k_raw, tot)             # 真实计数,精确
        m_lo, m_hi = analysis.wilson_ci(round(f_mit * tot), tot)  # 伪计数,名义
        rows.append({"label": r["label"], "pattern": bits,
                     "raw_fid": round(f_raw, 4),
                     "raw_fid_ci": [round(r_lo, 4), round(r_hi, 4)],
                     "mit_fid": round(f_mit, 4),
                     "mit_fid_ci_nominal": [round(m_lo, 4), round(m_hi, 4)],
                     "clipped": round(clip, 4)})
        print(f"  {r['label']:12s}: 恢复保真 {f_raw:.4f} → {f_mit:.4f} "
              f"(截断 {clip:.3f})")
    if rows:
        mean_raw = sum(r["raw_fid"] for r in rows) / len(rows)
        mean_mit = sum(r["mit_fid"] for r in rows) / len(rows)
        print(f"  均值: {mean_raw:.4f} → {mean_mit:.4f};残差 {1 - mean_mit:.4f} "
              f"= 上下文效应等模型偏差的实测上界")
    raw_pool = analysis.pool_ci([r["raw_fid"] for r in rows]) if rows else {"n": 0}
    mit_pool = analysis.pool_ci([r["mit_fid"] for r in rows]) if rows else {"n": 0}
    prov = provenance.stamp(
        anchor=a, qubits=reg, tids=provenance.tids_of(recs), chip=CHIP,
        extra={"selection": "anchmod.latest chains['7'] (当日最优 d7 链)",
               "cal_source": cal_src, "n_patterns": len(patterns),
               "shots_per_task": SHOTS,
               "ci": "raw_fid=Wilson 95% 精确;mit_fid_ci_nominal=伪计数 shot-noise "
                     "名义下界(反演放大真实方差,亦不含上下文效应)"})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"validation_{RUN}.json").write_text(
        json.dumps({"provenance": prov,
                    "aggregate": {"raw_fid": raw_pool, "mit_fid": mit_pool},
                    "rows": rows}, ensure_ascii=False, indent=2))
    print(f"已写 {OUT / f'validation_{RUN}.json'}")


if __name__ == "__main__":
    main()
