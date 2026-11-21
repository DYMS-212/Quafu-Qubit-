"""e20 · 冻结锚点 vs 当日 η 锚点:同批 A/B(切开"机器漂移"与"选址指标好坏")。

问题:当日 η 选的链,在 repcode 上未必优于冻结参考链(07-04 实测:冻结链 d7L1=0.0044
反优于当日链 0.0105)。本实验同一天、同一批次**交错**跑冻结链与当日链的同一 repcode
(d=3/5/7 bit L1),直接量"当日选址增益 = 冻结err/当日err":
  - 增益 >1 稳定 → η 选址确实有用;
  - 增益 <1(冻结反而好)反复出现 → η 不是 repcode 的充分选择指标(e21 找更好的)。
两臂同批交错提交(避免日内漂移污染),附 Wilson CI,多天积累看方向是否稳定。
冻结锚点是手工冻结的漂移标尺(只有 chains、无 star,属正常)。

运行: ./qrun experiments/e20_anchor_ab.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, anchors as anchmod, client, repcode, provenance  # noqa: E402

CHIP = "Baihua"
RUN = time.strftime("%Y-%m-%d_%H%M")
EXP = f"e20_anchor_ab/{RUN}"
OUT = Path(__file__).resolve().parent.parent / "results" / "e20_anchor_ab"
SHOTS = 8192
DS = [3, 5, 7]


def main():
    print(f"后端状态: {client.status()}")
    frozen = anchmod.load(CHIP)          # 冻结锚点(漂移标尺)
    daily = anchmod.latest(CHIP)         # 当日 η 锚点
    arms = {"frozen": frozen, "daily": daily}

    specs = []
    for d in DS:                          # 同 d 的两臂相邻 → 同批交错
        for name, a in arms.items():
            chain = a["chains"][str(d)]
            specs.append(dict(circuit=repcode.build_t1(d, 1), chip=CHIP, shots=SHOTS,
                              name=f"e20_rep_d{d}_{name}", exp=EXP,
                              label=f"rep_d{d}L1_{name}",
                              options={"target_qubits": chain},
                              meta={"kind": "rep", "d": d, "arm": name, "qubits": chain}))
    recs = client.run_batch(specs, timeout=1200)

    rows = []
    for r in recs:
        c = client.counts(r)
        if not c:
            continue
        m = r["meta"]
        err = repcode.logical_error_majority(c, m["d"], 1)
        sh = analysis.shots_of(c)
        rows.append({"kind": "rep", "d": m["d"], "arm": m["arm"], "qubits": m["qubits"],
                     "tid": r.get("tid"), "err": round(err, 6),
                     "ci": [round(x, 6) for x in analysis.wilson_ci(round(err * sh), sh)],
                     "shots": sh})

    print("\n===== 冻结 vs 当日锚点:repcode 选址增益(同批 A/B)=====")
    for d in DS:
        fr = next((r for r in rows if r["d"] == d and r["arm"] == "frozen"), None)
        da = next((r for r in rows if r["d"] == d and r["arm"] == "daily"), None)
        if fr and da:
            gain = fr["err"] / (da["err"] or 1e-9)
            print(f"  d={d}: 冻结={fr['err']:.5f} 当日={da['err']:.5f} "
                  f"→ 当日选址增益 {gain:.2f}× ({'当日更好' if gain > 1 else '冻结反而更好'})")

    OUT.mkdir(parents=True, exist_ok=True)
    summary = {"provenance": provenance.stamp(
                   anchor=daily, chip=CHIP, tids=provenance.tids_of(recs),
                   qubits={"frozen_chains": frozen.get("chains"),
                           "daily_chains": daily.get("chains")},
                   extra={"frozen_cal_time": frozen.get("chip_calibration_time"),
                          "note": "冻结锚点=手工冻结的漂移标尺(无star属正常);本实验只比链选址"}),
               "rows": rows}
    (OUT / f"summary_{RUN}.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"已写 {OUT / f'summary_{RUN}.json'}")


if __name__ == "__main__":
    main()
