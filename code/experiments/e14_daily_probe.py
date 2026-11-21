"""e14 · 每日体检(固定锚点,天天跑,漂移追踪 + 置信区间的数据地基)。

固定套件(13 个任务,位置永不变,保证跨日可比):
  1. 全芯片自校准 all0/all1(available 比特)         × 4096 shots
  2. ZZ syndrome 锚点 B0: data[20,22] anc[21],4 注错 × 2048
  3. repcode d=3 链 [8,7,6,5,4]: L0 / L1 / L0+中位注错 × 4096
  4. Bell 边 (67,68) 与 (73,74): ZZ / XX 两基          × 2048

产出:results/e14_probe/<日期>/*.json(原始)+ probe_log.jsonl 追加一行当日汇总。
运行: ./qrun experiments/e14_daily_probe.py [chip]   # 默认 Baihua
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, chip as chipmod, client, decoder, provenance, repcode  # noqa: E402
from qfab.qasm import QASM, bell_pair  # noqa: E402

CHIP = sys.argv[1] if len(sys.argv) > 1 else "Baihua"
MODE = sys.argv[2].lower() if len(sys.argv) > 2 else "full"  # full=每日全套体检 / light=贴 headline 每轮的轻量漂移快照(仅自校准)
TODAY = time.strftime("%Y-%m-%d")
EXP = f"e14_probe/{TODAY}"
OUT = Path(__file__).resolve().parent.parent / "results" / "e14_probe"

# ---- 固定锚点(勿改,跨日可比性所系) ----
SYN_DATA, SYN_ANC = [20, 22], 21          # B0 syndrome 锚点
REP_CHAIN = [8, 7, 6, 5, 4]               # d=3: data[8,6,4] anc[7,5]
BELL_EDGES = [(67, 68), (73, 74)]         # 干净直线 67-74 的两端


def syndrome_circ(case: str) -> str:
    qc = QASM(3)
    if "X0" in case:
        qc.x(0)
    if "X1" in case:
        qc.x(1)
    qc.barrier()
    qc.cx(0, 2).cx(1, 2)
    qc.barrier().measure_all()
    return qc.qasm()


def build_specs(av: list[int], mode: str = "full") -> list[dict]:
    specs = []
    n = len(av)
    # 1. 自校准(full/light 都跑——这是漂移序列的核心)
    qc0 = QASM(n); qc0.measure_all()
    qc1 = QASM(n)
    for i in range(n):
        qc1.x(i)
    qc1.measure_all()
    for tag, qc in (("prep0", qc0), ("prep1", qc1)):
        specs.append(dict(circuit=qc.qasm(), chip=CHIP, shots=4096,
                          name=f"probe_selfcal_{tag}", exp=EXP,
                          label=f"selfcal_{tag}",
                          options={"target_qubits": av},
                          meta={"kind": "selfcal", "prep": tag, "qubits": av}))
    if mode == "light":
        return specs   # 轻量:只自校准 2 任务,贴 headline 每轮,产出日内读出漂移序列
    # 2. syndrome 锚点
    for case in ("none", "X0", "X1", "X0X1"):
        specs.append(dict(circuit=syndrome_circ(case), chip=CHIP, shots=2048,
                          name=f"probe_syn_{case}", exp=EXP, label=f"syn_{case}",
                          options={"target_qubits": [SYN_DATA[0], SYN_DATA[1], SYN_ANC]},
                          meta={"kind": "syndrome", "case": case}))
    # 3. repcode d=3
    for tag, bit, inj in (("L0", 0, []), ("L1", 1, []), ("L0_injmid", 0, [1])):
        specs.append(dict(circuit=repcode.build_t1(3, bit, inj), chip=CHIP,
                          shots=4096, name=f"probe_rep_{tag}", exp=EXP,
                          label=f"rep_{tag}",
                          options={"target_qubits": REP_CHAIN},
                          meta={"kind": "repcode", "d": 3, "bit": bit, "inj": inj}))
    # 4. Bell 边
    for (a, b) in BELL_EDGES:
        for basis in ("ZZ", "XX"):
            specs.append(dict(circuit=bell_pair(basis).qasm(), chip=CHIP,
                              shots=2048, name=f"probe_bell_{a}_{b}_{basis}",
                              exp=EXP, label=f"bell_{a}_{b}_{basis}",
                              options={"target_qubits": [a, b]},
                              meta={"kind": "bell", "edge": [a, b], "basis": basis}))
    return specs


def summarize(recs: list[dict], av: list[int], cal_time) -> dict:
    by = {r["label"]: r for r in recs}

    def cnt(label):
        return client.counts(by[label]) if label in by else {}

    day = {"date": TODAY, "time": time.strftime("%H:%M"),  # 同日多跑=日内漂移数据
           "chip": CHIP, "vendor_calibration_time": cal_time,
           "n_available": len(av)}
    # 自校准
    c0, c1 = cnt("selfcal_prep0"), cnt("selfcal_prep1")
    if c0 and c1:
        n = len(av)
        asg = analysis.assignment_from_p1(analysis.bit_p1(c0, n), analysis.bit_p1(c1, n))
        etas = sorted(a["eta"] for a in asg)
        day["selfcal"] = {
            "eta_median": round(etas[len(etas) // 2], 4),
            "n_eta_below_0.5": sum(e < 0.5 for e in etas),
            "anchor_qubits": {str(q): {k: round(asg[av.index(q)][k], 4)
                                       for k in ("F0", "F1", "eta")}
                              for q in [*SYN_DATA, SYN_ANC, *REP_CHAIN,
                                        *[x for e in BELL_EDGES for x in e]]
                              if q in av},
        }
    # syndrome:P(anc=期望)(+ Wilson CI,计数直接取自原始 counts)
    expect = {"none": 0, "X0": 1, "X1": 1, "X0X1": 0}
    syn, syn_ci = {}, {}
    for case, exp_s in expect.items():
        c = cnt(f"syn_{case}")
        if c:
            n = analysis.shots_of(c)
            p1 = analysis.bit_p1(c, 3)[2]
            rate = p1 if exp_s else 1 - p1
            syn[case] = round(rate, 4)
            lo, hi = analysis.wilson_ci(round(rate * n), n)
            syn_ci[case] = [round(lo, 4), round(hi, 4)]
    if syn:
        day["syndrome_success"] = syn
        day["syndrome_success_ci"] = syn_ci
    # repcode(majority_err + Wilson CI)
    rep = {}
    for tag, bit in (("L0", 0), ("L1", 1), ("L0_injmid", 0)):
        c = cnt(f"rep_{tag}")
        if c:
            n = analysis.shots_of(c)
            me = repcode.logical_error_majority(c, 3, bit)
            lo, hi = analysis.wilson_ci(round(me * n), n)
            rep[tag] = {
                "majority_err": round(me, 5),
                "majority_err_ci": [round(lo, 5), round(hi, 5)],
                "detector_rates": [round(x, 4) for x in decoder.detector_rates(c, 3)],
            }
    if rep:
        day["repcode_d3"] = rep
    # bell(parity + CI:Wilson 建在 even-fraction 上再映射回 parity 尺度)
    bell, bell_ci = {}, {}
    for (a, b) in BELL_EDGES:
        for basis in ("ZZ", "XX"):
            c = cnt(f"bell_{a}_{b}_{basis}")
            if c:
                n = analysis.shots_of(c)
                even = sum(v for k, v in c.items()
                           if k.zfill(2).count("1") % 2 == 0)
                key = f"{a}-{b}_{basis}"
                bell[key] = round(2 * even / n - 1, 4)
                lo, hi = analysis.wilson_ci(even, n)
                bell_ci[key] = [round(2 * lo - 1, 4), round(2 * hi - 1, 4)]
    if bell:
        day["bell_parity"] = bell
        day["bell_parity_ci"] = bell_ci
    # 溯源戳:冻结硬编码锚点、13 个任务 id、完整可用比特(自校准覆盖面)。
    day["provenance"] = provenance.stamp(
        anchor=None, qubits=av, tids=provenance.tids_of(recs), chip=CHIP,
        extra={"anchor": "frozen_hardcoded",
               "chip_calibration_time": cal_time,
               "probe_anchor_qubits": {
                   "syndrome": [SYN_DATA[0], SYN_DATA[1], SYN_ANC],
                   "rep_chain_d3": REP_CHAIN,
                   "bell_edges": [list(e) for e in BELL_EDGES]}})
    return day


def main():
    st = client.status()
    print(f"后端状态: {st}")
    if not isinstance(st.get(CHIP), int):
        print(f"⚠️ {CHIP} 当前不可用({st.get(CHIP)}),中止。")
        return
    ch = chipmod.load(CHIP, refresh=True)
    print(f"芯片: {ch.summary()}")
    av = ch.available_qubits()
    specs = build_specs(av, MODE)
    recs = client.run_batch(specs, timeout=1200)
    day = summarize(recs, av, ch.calibration_time)
    day["mode"] = MODE   # full=每日全套体检 / light=贴 headline 轮的漂移快照
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "probe_log.jsonl").open("a", encoding="utf-8") as f:
        f.write(json.dumps(day, ensure_ascii=False) + "\n")
    if MODE == "light":
        sc = day.get("selfcal", {})
        print(f"[light 漂移快照 {day.get('time')}] η中位={sc.get('eta_median')} "
              f"η<0.5={sc.get('n_eta_below_0.5')}")
    else:
        print("\n===== 当日体检汇总 =====")
        print(json.dumps(day, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
