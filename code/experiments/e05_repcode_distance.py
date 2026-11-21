"""e05 · repetition 码距离扫描(d=3..13 × bit/phase × L0/L1,带 fresh-ancilla syndrome)。

d 扫到锚点长链允许的最大值(d=13 → 25 比特):p_L(d) 曲线从"三个点"变成
"压制到地板/失效点的完整曲线"——距离扩展在哪儿失效、为什么,本身就是结果。

产出三类下游:
  - 距离压制曲线 p_L(d)(多天积累误差棒)
  - e06 检错后选(离线吃本数据)
  - e12 解码器阶梯(离线吃本数据 + e14 自测权重)

运行: ./qrun experiments/e05_repcode_distance.py [chip]   # 默认 Baihua
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, anchors as anchmod, client, repcode  # noqa: E402

CHIP = sys.argv[1] if len(sys.argv) > 1 else "Baihua"
TODAY = time.strftime("%Y-%m-%d")
RUN = time.strftime("%Y-%m-%d_%H%M")   # 同日多轮互不覆盖
TAG = RUN if CHIP == "Baihua" else f"{CHIP}_{RUN}"  # 非主线机文件名带机名防串味
EXP = f"e05_repcode/{TAG}"
OUT = Path(__file__).resolve().parent.parent / "results" / "e05_repcode"
SHOTS = 8192

_A = anchmod.latest(CHIP)
CHAINS = {int(k): v for k, v in _A["chains"].items()}


def main():
    print(f"后端状态: {client.status()}")
    specs = []
    for d, chain in CHAINS.items():
        for phase in (False, True):
            for bit in (0, 1):
                code = "phase" if phase else "bit"
                specs.append(dict(
                    circuit=repcode.build_t1(d, bit, phase=phase),
                    chip=CHIP, shots=SHOTS,
                    name=f"rep_{code}_d{d}L{bit}", exp=EXP,
                    label=f"{code}_d{d}_L{bit}",
                    options={"target_qubits": chain},
                    meta={"d": d, "bit": bit, "phase": phase, "code": code}))
    print(f"共 {len(specs)} 条线路 × {SHOTS} shots")
    recs = client.run_batch(specs, timeout=1200)

    rows = []
    print("\n===== p_L(d) 多数表决(检错后选/解码器见 e06/e12 离线分析)=====")
    for r in recs:
        c = client.counts(r)
        if not c:
            continue
        m = r["meta"]
        err = repcode.logical_error_majority(c, m["d"], m["bit"])
        k = round(err * analysis.shots_of(c))
        ci = analysis.wilson_ci(k, analysis.shots_of(c))
        rows.append({**m, "chip": CHIP, "chain": CHAINS.get(m["d"]), "tid": r.get("tid"),
                     "anchor_date": _A.get("date"),
                     "cal_time": _A.get("chip_calibration_time"),
                     "anchor_source": _A.get("generated_from"),
                     "err": round(err, 6), "ci": [round(x, 6) for x in ci],
                     "shots": analysis.shots_of(c)})
        print(f"  {m['code']:5s} d={m['d']} L{m['bit']}: {err:.5f} "
              f"[{ci[0]:.5f}, {ci[1]:.5f}]")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"summary_{TAG}.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2))
    print(f"已写 {OUT / f'summary_{TAG}.json'}")


if __name__ == "__main__":
    main()
