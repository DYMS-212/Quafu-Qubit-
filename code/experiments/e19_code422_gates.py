"""e19 · [[4,2,2]] 门级对决:逻辑门深度扫描(Gottesman 式编码 vs 裸电路判据)。

[[4,2,2]] 的 transversal H(H⊗4)= 双逻辑 H · 逻辑 SWAP(非平凡逻辑门)。
序列:|00⟩_L → (H̄)^k → 测量:
  k 偶 → 末态 |00⟩_L,Z 基读 Z̄₁=b0⊕b2, Z̄₂=b0⊕b1
  k 奇 → 末态 |++⟩_L,X 基读(测量前再加 H⊗4)X̄₁=b0⊕b1, X̄₂=b0⊕b2
  每层之间 barrier(防编译器把 H 对消——否则白测)。
后选 stabilizer 总奇偶。对照:2 个最好物理比特打同样 k 层 H(同样 barrier)。
判据:per-logical-qubit 后选错误率随 k 的斜率 vs 裸比特斜率。

**主指标 = 每层错误率 r**:多日数据证明"交叉点"由两臂截距
(制备/读出偏置)决定、随裸基线日间波动(0.02-0.09)漂移,不稳健;斜率
(RB 式衰减拟合 ln(1-2p) ~ k,见 analysis.fit_decay_per_layer)只反映每层门
质量,跨日可比。交叉点降级为条件性结果(裸截距差的日子才出现)。

选址:当日最优星形(anchors.daily)。
运行: ./qrun experiments/e19_code422_gates.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, anchors as anchmod, client, provenance  # noqa: E402
from qfab.qasm import QASM  # noqa: E402

CHIP = "Baihua"
RUN = time.strftime("%Y-%m-%d_%H%M")
EXP = f"e19_gates/{RUN}"
OUT = Path(__file__).resolve().parent.parent / "results" / "e19_gates"
SHOTS = 4096
KS = [0, 1, 2, 4, 6, 8, 12, 16, 20]   # 斜率(主指标)拟合需要长杠杆;
                                       # 多日回填显示 k≤12 时编码臂斜率≈噪声水平(见 tools/e19_slopes.py)


def circ_encoded(k: int) -> str:
    qc = QASM(4)
    qc.h(0).cx(0, 1).cx(0, 2).cx(0, 3)     # |00>_L
    qc.barrier()
    for _ in range(k):                      # k 层 transversal H = 逻辑 (H̄H̄·SWAP)^k
        for q in range(4):
            qc.h(q)
        qc.barrier()
    if k % 2 == 1:                          # 末态 |++>_L → X 基读出
        for q in range(4):
            qc.h(q)
    qc.barrier().measure_all()
    return qc.qasm()


def circ_bare(k: int) -> str:
    qc = QASM(2)
    qc.barrier()
    for _ in range(k):
        qc.h(0).h(1)
        qc.barrier()
    if k % 2 == 1:
        qc.h(0).h(1)                        # |+> → 读回 Z 基
    qc.barrier().measure_all()
    return qc.qasm()


def eval_encoded(count: dict, k: int) -> dict:
    """后选总奇偶;逻辑读出(全部期望 0):k偶 Z̄=(b0⊕b2, b0⊕b1),k奇 X̄=(b0⊕b1, b0⊕b2)。"""
    tot = acc = e = 0
    for key, v in count.items():
        b = analysis.bits_of(key, 4)
        tot += v
        if (b[0] ^ b[1] ^ b[2] ^ b[3]) == 0:
            acc += v
            if k % 2 == 0:
                l1, l2 = b[0] ^ b[2], b[0] ^ b[1]
            else:
                l1, l2 = b[0] ^ b[1], b[0] ^ b[2]
            e += v * ((l1 != 0) + (l2 != 0))
    return {"acceptance": acc / (tot or 1),
            "acceptance_ci": analysis.wilson_ci(acc, tot or 1),
            "ps_err_perq": e / (2 * (acc or 1)),
            "ps_err_perq_ci": analysis.wilson_ci(e, 2 * (acc or 1))}


def main():
    print(f"后端状态: {client.status()}")
    a = anchmod.latest(CHIP)
    star = a["star"]
    tq_enc = [star["center"], *star["leaves"]]
    bare = a["phys_baseline"]
    print(f"星形: {tq_enc},裸对照: {bare}")

    ss = tq_enc[:2]  # 同substrate裸对照:星形的中心+首叶(与编码同一批物理比特)
    print(f"同substrate裸对照(星形子集): {ss}")
    specs = []
    for k in KS:                            # enc / bare(最强物理) / bare_ss(同substrate)交错
        specs.append(dict(circuit=circ_encoded(k), chip=CHIP, shots=SHOTS,
                          name=f"g422_enc_k{k}", exp=EXP, label=f"enc_k{k}",
                          options={"target_qubits": tq_enc},
                          meta={"kind": "enc", "k": k, "star": tq_enc}))
        specs.append(dict(circuit=circ_bare(k), chip=CHIP, shots=SHOTS,
                          name=f"g422_bare_k{k}", exp=EXP, label=f"bare_k{k}",
                          options={"target_qubits": bare},
                          meta={"kind": "bare", "k": k, "qubits": bare}))
        specs.append(dict(circuit=circ_bare(k), chip=CHIP, shots=SHOTS,
                          name=f"g422_baress_k{k}", exp=EXP, label=f"baress_k{k}",
                          options={"target_qubits": ss},
                          meta={"kind": "bare_ss", "k": k, "qubits": ss}))
    recs = client.run_batch(specs, timeout=1200)

    rows = []
    print("\n===== 逻辑门深度扫描(per-qubit 错误率,全部期望读 0)=====")
    print(f"{'k':>3s} {'编码(后选)':>12s} {'接受率':>8s} {'裸比特':>10s}")
    def bare_metrics(rec):
        c = client.counts(rec)
        p1 = analysis.bit_p1(c, 2)
        tot_b = analysis.shots_of(c)
        ones = int(round((p1[0] + p1[1]) * tot_b))  # 两裸比特错误合计计数
        return (p1[0] + p1[1]) / 2, analysis.wilson_ci(ones, 2 * (tot_b or 1))

    for k in KS:
        enc = next((r for r in recs if r["meta"]["kind"] == "enc"
                    and r["meta"]["k"] == k and client.counts(r)), None)
        row = {"k": k}
        if enc:
            row.update(eval_encoded(client.counts(enc), k))
        for tag in ("bare", "bare_ss"):     # bare=最强物理基线;bare_ss=同substrate
            rec = next((r for r in recs if r["meta"]["kind"] == tag
                        and r["meta"]["k"] == k and client.counts(r)), None)
            if rec:
                err, ci = bare_metrics(rec)
                row[f"{tag}_err_perq"] = err
                row[f"{tag}_err_perq_ci"] = ci
        rows.append({kk: (round(v, 5) if isinstance(v, float)
                          else [round(x, 5) for x in v] if isinstance(v, (tuple, list))
                          else v)
                     for kk, v in row.items()})
        print(f"{k:>3d} {row.get('ps_err_perq', float('nan')):>12.5f} "
              f"{row.get('acceptance', float('nan')):>8.3f} "
              f"{row.get('bare_err_perq', float('nan')):>10.5f}")
    # 每层错误率(主指标):RB 式衰减拟合,截距漂移不影响
    slopes = {}
    for key, field in (("enc", "ps_err_perq"), ("bare", "bare_err_perq"),
                       ("bare_ss", "bare_ss_err_perq")):
        fit = analysis.fit_decay_per_layer(
            [(r["k"], r.get(field)) for r in rows])
        slopes[key] = {kk: (round(v, 6) if isinstance(v, float) else v)
                       for kk, v in fit.items()}
    if slopes["enc"].get("r") and slopes["bare"].get("r"):
        slopes["ratio_bare_over_enc"] = round(
            slopes["bare"]["r"] / max(slopes["enc"]["r"], 1e-9), 3)
    print("\n每层错误率 r(拟合 ln(1-2p)~k;主指标,跨日可比):")
    for key in ("enc", "bare", "bare_ss"):
        f = slopes[key]
        if f.get("r") is not None:
            print(f"  {key:8s} r={f['r']:.5f} CI[{f['r_lo']:.5f},{f['r_hi']:.5f}] "
                  f"截距损耗≈{f['prep_offset']:.4f} (n={f['n']})")
    if "ratio_bare_over_enc" in slopes:
        print(f"  斜率比 bare/enc = {slopes['ratio_bare_over_enc']}×")

    OUT.mkdir(parents=True, exist_ok=True)
    summary = {"provenance": provenance.stamp(
                   anchor=a, qubits={"star": tq_enc, "bare_best": bare, "bare_ss": ss},
                   tids=provenance.tids_of(recs), chip=CHIP),
               "slopes": slopes,
               "rows": rows}
    (OUT / f"summary_{RUN}.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"已写 {OUT / f'summary_{RUN}.json'}")


if __name__ == "__main__":
    main()
