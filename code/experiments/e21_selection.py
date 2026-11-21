"""e21 · 选址指标 bake-off:哪个指标最能预测 d7 repcode 性能?

背景:vendor gate_fid 近乎均匀(std 0.0045,无区分度);若选址只看读出 η、不看 CZ
门质量,可能选到"读出好但门差"的链(编码/深门电路吃亏)。本实验用同一份当日自测
η + 芯片 CZ,按不同指标各选一条 d7 链,同批跑 repcode d7 bit L1,比谁的链错误率最低:
  - readout   : 只看读出 η(当前默认)
  - composite : η(节点)× CZ(沿途边)—— 兼顾读出与门
  - frozen    : 冻结参考链(基线)
胜者(多天稳定最低)= 应接进 anchors 的选址指标。若 readout 与 composite 选出同一条链,
即"η 已隐含选到 CZ 好的链",也是一条结论。(自测 CZ 代理版见 e03 Bell parity,后续接入。)

运行: ./qrun experiments/e21_selection.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, anchors as anchmod, chip as chipmod, client, repcode, provenance  # noqa: E402

CHIP = "Baihua"
RUN = time.strftime("%Y-%m-%d_%H%M")
EXP = f"e21_selection/{RUN}"
PROJ = Path(__file__).resolve().parent.parent
OUT = PROJ / "results" / "e21_selection"
SHOTS = 8192
MODES = ["readout", "composite"]


def load_today_eta() -> dict[int, float]:
    d = PROJ / "results" / "e14_probe" / time.strftime("%Y-%m-%d")
    rec0 = rec1 = None
    for f in d.glob("selfcal_prep0__*.json"):
        rec0 = json.loads(f.read_text())
    for f in d.glob("selfcal_prep1__*.json"):
        rec1 = json.loads(f.read_text())
    if not (rec0 and rec1):
        raise RuntimeError("缺当日自校准,先跑 e14_daily_probe")
    qubits = rec0["meta"]["qubits"]
    asg = analysis.assignment_from_p1(
        analysis.bit_p1(rec0["result"]["count"], len(qubits)),
        analysis.bit_p1(rec1["result"]["count"], len(qubits)))
    return {q: asg[i]["eta"] for i, q in enumerate(qubits)}


def main():
    print(f"后端状态: {client.status()}")
    ch = chipmod.load(CHIP, refresh=True)
    eta = load_today_eta()
    chains = {mode: anchmod._generate(ch, eta, mode)["chains"].get("7") for mode in MODES}
    chains["frozen"] = anchmod.load(CHIP)["chains"].get("7")

    # 相同链只跑一次,记录哪些指标共享它
    uniq: dict = {}
    for label, chain in chains.items():
        uniq.setdefault(tuple(chain) if chain else None, []).append(label)
    print("各指标选出的 d7 链:")
    for key, labels in uniq.items():
        print(f"  {'+'.join(labels)}: {list(key) if key else None}")

    specs = []
    for key, labels in uniq.items():
        if not key:
            continue
        tag = "+".join(labels)
        specs.append(dict(circuit=repcode.build_t1(7, 1), chip=CHIP, shots=SHOTS,
                          name=f"e21_{tag}", exp=EXP, label=f"d7L1_{tag}",
                          options={"target_qubits": list(key)},
                          meta={"metrics": labels, "chain": list(key)}))
    recs = client.run_batch(specs, timeout=1200)

    rows = []
    for r in recs:
        c = client.counts(r)
        if not c:
            continue
        m = r["meta"]
        err = repcode.logical_error_majority(c, 7, 1)
        sh = analysis.shots_of(c)
        rows.append({"metrics": m["metrics"], "chain": m["chain"], "tid": r.get("tid"),
                     "err": round(err, 6),
                     "ci": [round(x, 6) for x in analysis.wilson_ci(round(err * sh), sh)],
                     "shots": sh})

    rows.sort(key=lambda r: r["err"])
    print("\n===== 选址指标 bake-off(d7 repcode bit L1,越低越好)=====")
    for r in rows:
        print(f"  {'+'.join(r['metrics']):22s} err={r['err']:.5f} "
              f"[{r['ci'][0]:.5f},{r['ci'][1]:.5f}]")
    if rows:
        print(f"  → 今日最佳指标: {'+'.join(rows[0]['metrics'])}")

    OUT.mkdir(parents=True, exist_ok=True)
    summary = {"provenance": provenance.stamp(
                   anchor=None, chip=CHIP, tids=provenance.tids_of(recs),
                   qubits=chains,
                   extra={"modes": MODES, "chip_calibration_time": ch.calibration_time,
                          "measured_date": time.strftime("%Y-%m-%d")}),
               "rows": rows}
    (OUT / f"summary_{RUN}.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"已写 {OUT / f'summary_{RUN}.json'}")


if __name__ == "__main__":
    main()
