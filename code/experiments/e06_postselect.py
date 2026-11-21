"""e06 · 检错后选:接受率-错误率权衡(离线,吃 e05 数据,零额度)。

对每条 e05 记录:syndrome 全零后选 → 后选逻辑错误率 + 接受率,对照 raw 多数表决。
铁律:后选错误率必须与接受率一并报告。

运行: ./qrun experiments/e06_postselect.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, provenance, repcode  # noqa: E402

PROJ = Path(__file__).resolve().parent.parent
OUT = PROJ / "results" / "e06_postselect"


def postselect(count: dict, d: int, bit: int) -> dict:
    tot = acc = err = 0
    for key, v in count.items():
        data, synd = repcode.split_bits(key, d)
        tot += v
        if not any(synd):
            acc += v
            if repcode.majority(data) != bit:
                err += v
    return {"acceptance": acc / (tot or 1), "ps_err": err / (acc or 1),
            "tot": tot, "acc": acc, "err": err}


def main():
    days = sorted(d for d in (PROJ / "results" / "e05_repcode").iterdir()
                  if d.is_dir() and d.name[:2] == "20")  # 只吃 Baihua 主线(日期名目录)
    day = days[-1]
    print(f"数据取自 e05 {day.name}")
    rows = []
    used_recs = []
    print(f"{'code':6s} {'d':>2s} {'L':>2s} {'raw_maj':>9s} {'后选err':>9s} {'接受率':>7s}")
    for f in sorted(day.glob("*_d*_L*__*.json")):
        rec = json.loads(f.read_text())
        c = (rec.get("result") or {}).get("count")
        if not c:
            continue
        used_recs.append(rec)
        m = rec["meta"]
        chain = rec["task"]["options"]["target_qubits"]
        raw = repcode.logical_error_majority(c, m["d"], m["bit"])
        ps = postselect(c, m["d"], m["bit"])
        tot, acc, err = ps["tot"], ps["acc"], ps["err"]
        raw_lo, raw_hi = analysis.wilson_ci(round(raw * tot), tot)
        acc_lo, acc_hi = analysis.wilson_ci(acc, tot)          # 接受率 = acc/shots
        pse_lo, pse_hi = analysis.wilson_ci(err, acc)          # 后选err = err/acc
        rows.append({**m, "chain": chain,
                     "raw_majority_err": round(raw, 6),
                     "raw_majority_err_ci": [round(raw_lo, 6), round(raw_hi, 6)],
                     "ps_err": round(ps["ps_err"], 6),
                     "ps_err_ci": [round(pse_lo, 6), round(pse_hi, 6)],
                     "acceptance": round(ps["acceptance"], 4),
                     "acceptance_ci": [round(acc_lo, 4), round(acc_hi, 4)],
                     "shots": tot})
        print(f"{m['code']:6s} {m['d']:>2d} {m['bit']:>2d} {raw:>9.5f} "
              f"{ps['ps_err']:>9.5f} {ps['acceptance']:>7.3f}")
    chip = used_recs[0].get("chip", "Baihua") if used_recs else "Baihua"
    prov = provenance.stamp(
        chip=chip, tids=provenance.tids_of(used_recs),
        extra={"kind": "offline_reanalysis",
               "source_exp": "e05_repcode", "source_dir": day.name,
               "source_path": str(day), "n_records": len(used_recs),
               "note": "每行物理链见 row['chain']",
               "ci": "Wilson 95%:acceptance=acc/shots, ps_err=err/acc, raw=err/shots"})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"summary_{day.name}.json").write_text(
        json.dumps({"provenance": prov, "rows": rows}, ensure_ascii=False, indent=2))
    print(f"已写 {OUT / f'summary_{day.name}.json'}")


if __name__ == "__main__":
    main()
