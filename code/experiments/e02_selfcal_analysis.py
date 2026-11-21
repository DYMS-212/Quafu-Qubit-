"""e02 · 全芯片读出自校准分析 + 跨日漂移(离线,吃 e14 体检数据)。

产出:
  - 最新一天全芯片 assignment matrix(F0/F1/η/bias)+ 厂家发布值并列对照:
    厂家发布其标定快照、我们发布当日自测,两边各记时间戳,差异归因于"快照 vs 实时/上下文"(不作评判)
  - 跨日漂移时间线(η 中位、差比特数、锚点比特轨迹)

运行: ./qrun experiments/e02_selfcal_analysis.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, anchors  # noqa: E402

PROJ = Path(__file__).resolve().parent.parent
OUT = PROJ / "results" / "e02_selfcal"


def day_assignment(day: Path):
    """当天体检的读出 assignment(跳过队列超时未补收的记录,详见 anchors)。"""
    return anchors.assignment_from_probe_dir(day)


def vendor_readout(chip: str, date: str):
    """厂家发布的读出标定快照 + 其标定时刻(他们测的、他们发布的值)。"""
    snap = PROJ / "results" / "chip_snapshots" / f"chip_{chip}_{date}.json"
    if not snap.exists():
        return {}, None
    info = json.loads(snap.read_text())
    out = {}
    for _, qd in info.get("qubits_info", {}).items():
        vals = [v for k, v in qd.items()
                if "readout" in k and "fidelity" in k and isinstance(v, (int, float))]
        if vals:
            out[qd["index"]] = sum(vals) / len(vals)
    return out, info.get("calibration_time")


def main():
    days = sorted(d for d in (PROJ / "results" / "e14_probe").iterdir() if d.is_dir())
    OUT.mkdir(parents=True, exist_ok=True)

    print("===== 跨日漂移时间线 =====")
    timeline = []
    for day in days:
        qubits, asg = day_assignment(day)
        if not asg:
            continue
        etas = sorted(a["eta"] for a in asg)
        row = {"date": day.name,
               "eta_median": round(etas[len(etas) // 2], 4),
               "n_eta_below_0.5": sum(e < 0.5 for e in etas),
               "n_qubits": len(etas)}
        timeline.append(row)
        print(f"  {row['date']}: η中位={row['eta_median']}, "
              f"差比特(η<0.5)={row['n_eta_below_0.5']}/{row['n_qubits']}")

    # 最新有效一天:完整 assignment + vendor 对照(快照-实测差距)
    qubits = asg = None
    for day in reversed(days):
        qubits, asg = day_assignment(day)
        if asg:
            days = days[:days.index(day) + 1]
            break
    if not asg:
        print("没有任何有效自校准数据,先跑 e14")
        return
    ven, ven_caltime = vendor_readout("Baihua", days[-1].name)
    full = []
    for i, q in enumerate(qubits):
        row = {"qubit": q, **{k: round(v, 4) for k, v in asg[i].items()},
               "vendor": round(ven.get(q, float("nan")), 4)}
        row["vendor_minus_meas"] = round(row["vendor"] - row["F_mean"], 4) \
            if q in ven else None
        full.append(row)
    GAP_TH = 0.15  # 仅"差异最大一批"的展示阈值,非合格线/评判线
    # 厂家发布了有效值(>0)才谈"数值差异";发布 0/未标定的单列为"覆盖差异",不混为一谈
    valued = [r for r in full if (r["vendor"] or 0) > 0 and r["vendor_minus_meas"] is not None]
    gaps = sorted((r for r in valued if abs(r["vendor_minus_meas"]) > GAP_TH),
                  key=lambda r: -abs(r["vendor_minus_meas"]))
    ven_uncal = [r for r in full if not ((r["vendor"] or 0) > 0) and r["eta"] > 0.5]
    print(f"\n===== 厂家发布快照 vs 我们当日自测(差异=快照与实时/上下文的距离,非孰对孰错)=====")
    print(f"  厂家快照标定时刻={ven_caltime} | 我们自测日={days[-1].name}")
    print(f"  A) 双方都有值、|Δ|>{GAP_TH} 的比特,共 {len(gaps)} 个:")
    for r in gaps[:12]:
        print(f"     Q{r['qubit']}: 厂家={r['vendor']} 我们自测F_mean={r['F_mean']}"
              f" (F0={r['F0']} F1={r['F1']}) Δ={r['vendor_minus_meas']:+}")
    print(f"  B) 厂家快照标 0/未标定、我们自测却可用(η>0.5)的比特:共 {len(ven_uncal)} 个"
          f"(覆盖差异,非数值分歧)")

    (OUT / f"assignment_{days[-1].name}.json").write_text(
        json.dumps(full, ensure_ascii=False, indent=2))
    (OUT / "drift_timeline.json").write_text(
        json.dumps(timeline, ensure_ascii=False, indent=2))
    print(f"\n已写 assignment_{days[-1].name}.json + drift_timeline.json")

    # 每日自适应选址(headline 实验用;冻结锚点不受影响)
    eta_map = {q: asg[i]["eta"] for i, q in enumerate(qubits)}
    daily = anchors.daily_best("Baihua", eta_map, days[-1].name)
    print("\n===== 今日推荐锚点(自测 η 选址)=====")
    for k in ("chains", "syn_triple", "bell_edges", "load_edge",
              "phys_baseline", "star"):
        print(f"  {k}: {daily.get(k)}")


if __name__ == "__main__":
    main()
