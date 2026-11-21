"""e17 · 重复性/方差表征:误差棒的地基(评委问"你的误差棒哪来的"就用它答)。

同一线路在同一批次重复提交 R 次 → job-to-job 方差;
每个 job 内 shots 二项波动 → shot 噪声理论值;
两者对比:若 job 间方差 ≫ shot 噪声,说明存在批内漂移/编译抖动 →
   所有单次结果都必须打上这个"系统方差"的补丁(诚实统计的卖点)。

探针:冻结锚点 Bell 边 parity + repcode d3 L0(各重复 R 次)。
运行: ./qrun experiments/e17_variance.py [R]   # 默认 8
"""
from __future__ import annotations

import json
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import (analysis, anchors as anchmod, chip as chipmod,  # noqa: E402
                  client, provenance, repcode)
from qfab.qasm import bell_pair  # noqa: E402

CHIP = "Baihua"
RUN = time.strftime("%Y-%m-%d_%H%M")
EXP = f"e17_variance/{RUN}"
OUT = Path(__file__).resolve().parent.parent / "results" / "e17_variance"
SHOTS = 2048
R = int(sys.argv[1]) if len(sys.argv) > 1 else 8


def main():
    print(f"后端状态: {client.status()}")
    a = anchmod.load(CHIP)                     # 冻结锚点(跨日可比)
    edge = a["bell_edges"][0]
    chain = a["chains"]["3"]
    specs = []
    for r in range(R):
        specs.append(dict(circuit=bell_pair("ZZ").qasm(), chip=CHIP, shots=SHOTS,
                          name=f"var_bell_r{r}", exp=EXP, label=f"bell_r{r}",
                          options={"target_qubits": edge},
                          meta={"kind": "bell", "rep": r}))
        specs.append(dict(circuit=repcode.build_t1(3, 0), chip=CHIP, shots=SHOTS,
                          name=f"var_rep_r{r}", exp=EXP, label=f"rep_r{r}",
                          options={"target_qubits": chain},
                          meta={"kind": "rep", "rep": r}))
    recs = client.run_batch(specs, timeout=1200)

    vals: dict[str, list[float]] = {"bell": [], "rep": []}
    for r in recs:
        c = client.counts(r)
        if not c:
            continue
        if r["meta"]["kind"] == "bell":
            even = sum(v for k, v in c.items() if k.zfill(2).count("1") % 2 == 0)
            vals["bell"].append(even / analysis.shots_of(c))
        else:
            vals["rep"].append(repcode.logical_error_majority(c, 3, 0))
    out = {"run": RUN, "R": R, "shots": SHOTS, "edge": edge, "chain": chain}
    print("\n===== job-to-job 方差 vs shot 噪声 =====")
    for k, xs in vals.items():
        if len(xs) < 2:
            continue
        m = statistics.mean(xs)
        sd = statistics.stdev(xs)
        shot_sd = (m * (1 - m) / SHOTS) ** 0.5
        out[k] = {"mean": round(m, 5), "job_sd": round(sd, 5),
                  "shot_sd_theory": round(shot_sd, 5),
                  "excess_ratio": round(sd / shot_sd, 2) if shot_sd else None,
                  "pool_ci": analysis.pool_ci(xs),  # 跨 job 误差棒(mean/lo/hi/best/worst)
                  "values": [round(x, 5) for x in xs]}
        print(f"  {k}: mean={m:.5f} job间σ={sd:.5f} shot噪声σ={shot_sd:.5f} "
              f"→ 超额 {sd / shot_sd:.1f}×")
    # 溯源戳:冻结锚点身份(选址来源/锚点标定时刻)+ 本轮实机标定时刻 + 任务 id。
    try:
        run_cal = chipmod.load(CHIP, refresh=True).calibration_time
    except Exception:  # noqa: BLE001
        run_cal = None
    out["provenance"] = provenance.stamp(
        anchor=a, qubits=sorted(set(edge) | set(chain)),
        tids=provenance.tids_of(recs), chip=CHIP,
        extra={"anchor": "frozen", "site_selection": "frozen_anchor (anchors.load)",
               "run_calibration_time": run_cal})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"summary_{RUN}.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2))
    print(f"已写 {OUT / f'summary_{RUN}.json'}")


if __name__ == "__main__":
    main()
