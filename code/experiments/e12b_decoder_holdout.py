"""e12b · 解码器 holdout 验证(纯离线,配对 e12)。

回应质疑:e12 的 G_map_qfit 在同一份数据上拟合 q(物理翻转率)并评估,
可能高估。此处做严格 half-split:counts 按二项分布随机劈成 A/B 两半,
所有需要拟合的参数(q、ancilla 有效错误率)只用 A 半,错误率只在 B 半评估。
多次随机劈分取均值 ± 标准差 + 跨劈分 bootstrap 95% CI。

运行: ./qrun experiments/e12b_decoder_holdout.py
"""
from __future__ import annotations

import json
import random
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, anchors, decoder, provenance, repcode  # noqa: E402

PROJ = Path(__file__).resolve().parent.parent
OUT = PROJ / "results" / "e12_decoder"
N_SPLITS = 10
RUN = time.strftime("%Y-%m-%d_%H%M")


def split_half(count: dict, rng: random.Random) -> tuple[dict, dict]:
    a, b = {}, {}
    for k, v in count.items():
        na = sum(1 for _ in range(v) if rng.random() < 0.5)
        if na:
            a[k] = na
        if v - na:
            b[k] = v - na
    return a, b


def load_selfcal() -> tuple[dict[int, dict], str]:
    """最近一天有有效体检的自测 assignment(跳过队列超时未补收;今天全超时则回退最近有效日)。"""
    qubits, asg, dayname = anchors.latest_probe_assignment(PROJ / "results" / "e14_probe")
    if not qubits:
        raise RuntimeError("无任何有效体检自校准数据(先跑 e14_daily_probe + recollect)")
    return {q: asg[i] for i, q in enumerate(qubits)}, dayname


def main():
    cal, selfcal_dir = load_selfcal()
    results = []
    used_recs = []
    days = sorted(d for d in (PROJ / "results" / "e05_repcode").iterdir()
                  if d.is_dir() and d.name[:2] == "20")  # 只吃 Baihua 主线(日期名目录)
    print(f"数据取自 e05 {days[-1].name}")
    for f in sorted(days[-1].glob("*_d*_L*__*.json")):
        rec = json.loads(f.read_text())
        if not (rec.get("result") or {}).get("count"):
            continue
        used_recs.append(rec)
        d, bit = rec["meta"]["d"], rec["meta"]["bit"]
        chain = rec["task"]["options"]["target_qubits"]
        count = rec["result"]["count"]

        # 读出翻转先验(自测标定,独立数据,无需劈分)
        p_data, p_anc_cal = [], []
        for i, q in enumerate(chain):
            row = cal.get(q)
            if i % 2 == 0:
                p_data.append(max(1 - (row["F1"] if bit else row["F0"]), 1e-4))
            else:
                p_anc_cal.append(max(1 - row["F0"], 1e-4))

        maj_b, map_b = [], []
        rng = random.Random(7)
        for _ in range(N_SPLITS):
            half_a, half_b = split_half(count, rng)
            # 在 A 半拟合:ancilla 有效错误率 + 物理翻转率
            p_anc = decoder.p_anc_from_detectors(half_a, d, p_data)
            p1a = analysis.bit_p1(half_a, 2 * d - 1)
            q_fit = [max((p1a[2 * k] if bit == 0 else 1 - p1a[2 * k]) - p_data[k], 1e-4)
                     for k in range(d)]
            # 在 B 半评估
            maj_b.append(repcode.logical_error_majority(half_b, d, bit))
            map_b.append(decoder.logical_error_map(half_b, d, bit, q_fit, p_data, p_anc))

        maj_pool = analysis.pool_ci(maj_b)
        map_pool = analysis.pool_ci(map_b)
        row = {"code": rec["meta"].get("code", "bit"), "d": d, "bit": bit,
               "chain": chain, "n_splits": N_SPLITS,
               "majority": round(statistics.mean(maj_b), 5),
               "majority_sd": round(statistics.stdev(maj_b), 5),
               "majority_ci": [round(maj_pool["lo"], 5), round(maj_pool["hi"], 5)],
               "map_holdout": round(statistics.mean(map_b), 5),
               "map_holdout_sd": round(statistics.stdev(map_b), 5),
               "map_holdout_ci": [round(map_pool["lo"], 5), round(map_pool["hi"], 5)]}
        results.append(row)
        print(f"{row['code']} d={d} L{bit}: majority={row['majority']:.5f}±{row['majority_sd']:.5f}  "
              f"MAP(holdout)={row['map_holdout']:.5f}±{row['map_holdout_sd']:.5f}")

    chip = used_recs[0].get("chip", "Baihua") if used_recs else "Baihua"
    prov = provenance.stamp(
        chip=chip, tids=provenance.tids_of(used_recs),
        extra={"kind": "offline_holdout",
               "source_exp": "e05_repcode", "source_dir": days[-1].name,
               "source_path": str(days[-1]), "selfcal_dir": selfcal_dir,
               "n_records": len(used_recs), "n_splits": N_SPLITS,
               "note": "half-split:q 与 p_anc 只在 A 半拟合,错误率只在 B 半评估;"
                       "每行物理链见 row['chain'];_ci=跨劈分 bootstrap 95%"})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"holdout_{RUN}.json").write_text(
        json.dumps({"provenance": prov, "rows": results}, ensure_ascii=False, indent=2))
    print(f"\n已写 {OUT}/holdout_{RUN}.json")


if __name__ == "__main__":
    main()
