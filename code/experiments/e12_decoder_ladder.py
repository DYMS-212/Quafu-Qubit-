"""e12 · 解码器阶梯(纯离线,零真机成本)。

数据:e05 距离扫描(d=3/5/7 × bit/phase × L0/L1,带 fresh-ancilla syndrome,取最新一天)。
权重:e14 每日体检的全芯片自测 assignment matrix(取最新一天)。

解码器阶梯(逻辑错误率,越低越好):
  A majority        多数表决(忽略 syndrome)
  B soft-majority   校准加权表决(不用 syndrome)
  C MWPM-uniform    detector 图 MWPM,均匀权重
  D MWPM-selfcal    MWPM + 自测读出权重
  E MWPM-empirical  MWPM + 经验 detector 率交叉拟合(同码另一 prep,不循环)
  F/G map           精确边缘化 MAP(链 DP;F 零物理错先验,G 集内拟合)

运行: ./qrun experiments/e12_decoder_ladder.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, anchors, decoder, provenance, repcode  # noqa: E402

PROJ = Path(__file__).resolve().parent.parent
OUT = PROJ / "results" / "e12_decoder"

P_UNIFORM = 0.02
RUN = time.strftime("%Y-%m-%d_%H%M")


def load_selfcal() -> tuple[dict[int, dict], str]:
    """最近一天有有效体检的自测 assignment matrix(F0/F1/eta)。返回 (cal, 体检日目录名)。
    跳过队列超时未补收的记录;今天体检若全超时则自动回退到最近有效日(打印实际用的日期)。"""
    qubits, asg, dayname = anchors.latest_probe_assignment(PROJ / "results" / "e14_probe")
    if not qubits:
        raise RuntimeError("无任何有效体检自校准数据(先跑 e14_daily_probe + recollect)")
    print(f"自测权重取自体检 {dayname}")
    return {q: asg[i] for i, q in enumerate(qubits)}, dayname


def probs_for_chain(chain: list[int], prepared_bit: int, cal: dict[int, dict]):
    """state-conditional 读出翻转概率:data 期望 |bit>,ancilla 期望 |0>。"""
    d = (len(chain) + 1) // 2
    p_data, p_anc = [], []
    for i, q in enumerate(chain):
        row = cal.get(q)
        if row is None:  # 没自测数据就退回 uniform
            p = P_UNIFORM
        elif i % 2 == 0:  # data
            p = 1 - (row["F1"] if prepared_bit else row["F0"])
        else:  # ancilla 期望 0
            p = 1 - row["F0"]
        (p_data if i % 2 == 0 else p_anc).append(max(p, 1e-4))
    assert len(p_data) == d and len(p_anc) == d - 1
    return p_data, p_anc


def main():
    cal, selfcal_dir = load_selfcal()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    days = sorted(d for d in (PROJ / "results" / "e05_repcode").iterdir()
                  if d.is_dir() and d.name[:2] == "20")  # 只吃 Baihua 主线(日期名目录)
    files = sorted(days[-1].glob("*_d*_L*__*.json"))
    print(f"数据取自 e05 {days[-1].name}({len(files)} 条)")
    recs = [json.loads(f.read_text()) for f in files]
    recs = [r for r in recs if (r.get("result") or {}).get("count")]
    by_key = {(r["meta"]["code"], r["meta"]["d"], r["meta"]["bit"]): r for r in recs}
    for rec in recs:
        meta, count = rec["meta"], rec["result"]["count"]
        d, bit = meta["d"], meta["bit"]
        chain = rec["task"]["options"]["target_qubits"]
        shots = analysis.shots_of(count)

        p_data_cal, p_anc_cal = probs_for_chain(chain, bit, cal)
        p_data_u = [P_UNIFORM] * d
        p_anc_u = [P_UNIFORM] * (d - 1)

        # E:ancilla 有效错误率从"同码另一 prep"的 detector 率交叉拟合(不循环)
        other = by_key.get((meta["code"], d, 1 - bit))
        p_data_fit, _ = probs_for_chain(chain, 1 - bit, cal)
        p_anc_emp = decoder.p_anc_from_detectors(
            other["result"]["count"], d, p_data_fit) if other else p_anc_cal

        # F/G:精确边缘化 MAP。q=物理翻转率:F 用接近零的先验(完全不碰本集数据),
        # G 从本集边缘分布拟合(每比特 1 个标量,注意 in-sample,后续换 half-split 验证)
        q_floor = [1e-4] * d
        p1_obs = analysis.bit_p1(count, 2 * d - 1)
        q_fit = [max((p1_obs[2 * k] if bit == 0 else 1 - p1_obs[2 * k]) - p_data_cal[k], 1e-4)
                 for k in range(d)]

        decoders = {
            "A_majority": lambda c: repcode.logical_error_majority(c, d, bit),
            "B_softmaj": lambda c: decoder.logical_error_soft_majority(c, d, bit, p_data_cal),
            "C_mwpm_uniform": lambda c: decoder.logical_error(c, d, bit, p_data_u, p_anc_u),
            "D_mwpm_selfcal": lambda c: decoder.logical_error(c, d, bit, p_data_cal, p_anc_cal, soft=True),
            "E_mwpm_empirical": lambda c: decoder.logical_error(c, d, bit, p_data_cal, p_anc_emp, soft=True),
            "F_map_q0": lambda c: decoder.logical_error_map(c, d, bit, q_floor, p_data_cal, p_anc_emp),
            "G_map_qfit": lambda c: decoder.logical_error_map(c, d, bit, q_fit, p_data_cal, p_anc_emp),
        }
        row = {"code": meta["code"], "d": d, "bit": bit, "shots": shots, "chain": chain,
               "p_anc_selfcal": [round(p, 4) for p in p_anc_cal],
               "p_anc_empirical": [round(p, 4) for p in p_anc_emp]}
        for name, fn in decoders.items():
            ci = analysis.bootstrap_ci(count, fn, B=200, seed=42)
            row[name] = round(ci["value"], 6)
            row[name + "_ci"] = [round(ci["lo"], 6), round(ci["hi"], 6)]
        # G 的 q 在本集拟合又评估:bootstrap 冻结拟合参数 → CI 偏窄,真 holdout 见 e12b
        row["G_map_qfit_in_sample"] = True
        rows.append(row)
        print(f"{meta['code']} d={d} L{bit} shots={shots}")
        for name in decoders:
            lo, hi = row[name + "_ci"]
            print(f"    {name:16s} {row[name]:.5f}  [{lo:.5f}, {hi:.5f}]")

    chip = recs[0].get("chip", "Baihua") if recs else "Baihua"
    prov = provenance.stamp(
        chip=chip, tids=provenance.tids_of(recs),
        extra={"kind": "offline_reanalysis",
               "source_exp": "e05_repcode", "source_dir": days[-1].name,
               "source_path": str(days[-1]), "selfcal_dir": selfcal_dir,
               "n_records": len(recs), "p_uniform": P_UNIFORM,
               "note": "每行物理链见 row['chain'];各级 _ci=bootstrap 95%",
               "in_sample_ci_warning": "G_map_qfit 在本集拟合 q 又评估,bootstrap 冻结"
                                       "拟合参数 → 其 CI 偏窄(见 G_map_qfit_in_sample);"
                                       "真 holdout 验证见 e12b"})
    (OUT / f"summary_{RUN}.json").write_text(
        json.dumps({"provenance": prov, "rows": rows}, ensure_ascii=False, indent=2))
    print(f"\n已写 {OUT}/summary_{RUN}.json")

    # 汇总表
    print("\n===== 逻辑错误率对决(行=码,列=解码器)=====")
    hdr = (f"{'code':10s} {'A_maj':>9s} {'B_soft':>9s} {'C_mwpm_u':>9s} "
           f"{'D_mwpm_cal':>10s} {'E_mwpm_emp':>10s} {'F_map_q0':>9s} {'G_map_qfit':>10s}")
    print(hdr)
    for r in rows:
        print(f"{r['code'][:2]} d={r['d']} L{r['bit']:<3d} {r['A_majority']:>9.5f} {r['B_softmaj']:>9.5f} "
              f"{r['C_mwpm_uniform']:>9.5f} {r['D_mwpm_selfcal']:>10.5f} "
              f"{r['E_mwpm_empirical']:>10.5f} {r['F_map_q0']:>9.5f} {r['G_map_qfit']:>10.5f}")


if __name__ == "__main__":
    main()
