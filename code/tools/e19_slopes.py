"""e19 斜率回填:对全部历史 e19 summary 重算每层错误率 r(主指标)。

交叉点随裸基线截距日间漂移(0.02-0.09)不稳健;斜率只反映每层门质量。
本工具零额度,输出逐轮 r(enc/bare/bare_ss)与斜率比,并按日聚合 mean±CI,
供报告 §4.2 的多日统计与冻结判定使用。

运行: ./qrun tools/e19_slopes.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis  # noqa: E402

R = Path(__file__).resolve().parent.parent / "results" / "e19_gates"


def main():
    files = sorted(R.glob("summary_*.json")) if R.exists() else []
    if not files:
        print("无 e19 summary 可回填")
        return
    per_day: dict[str, dict[str, list]] = {}
    print(f"{'轮次':<22s} {'r_enc':>9s} {'r_bare':>9s} {'r_ss':>9s} {'bare/enc':>9s}")
    for f in files:
        data = json.loads(f.read_text())
        rows = data.get("rows", data) if isinstance(data, dict) else data
        day = f.name.replace("summary_", "")[:10]
        vals = {}
        for key, field in (("enc", "ps_err_perq"), ("bare", "bare_err_perq"),
                           ("bare_ss", "bare_ss_err_perq")):
            fit = analysis.fit_decay_per_layer(
                [(r["k"], r.get(field)) for r in rows])
            vals[key] = fit.get("r")
        ratio = (vals["bare"] / vals["enc"]
                 if vals.get("enc") and vals.get("bare") else None)
        fmt = lambda v: f"{v:9.5f}" if isinstance(v, float) else f"{'—':>9s}"  # noqa: E731
        print(f"{f.name.replace('summary_', '').replace('.json', ''):<22s} "
              f"{fmt(vals['enc'])} {fmt(vals['bare'])} {fmt(vals['bare_ss'])} "
              f"{fmt(ratio)}")
        d = per_day.setdefault(day, {"enc": [], "bare": [], "bare_ss": [], "ratio": []})
        for k in ("enc", "bare", "bare_ss"):
            if vals.get(k) is not None:
                d[k].append(vals[k])
        if ratio is not None:
            d["ratio"].append(ratio)

    print("\n按日聚合(mean±95%CI / n):")
    for day in sorted(per_day):
        parts = [day]
        for k in ("enc", "bare", "bare_ss", "ratio"):
            p = analysis.pool_ci(per_day[day][k])
            if p.get("n"):
                parts.append(f"{k}={p['mean']:.5f}[{p['lo']:.5f},{p['hi']:.5f}]n{p['n']}")
        print("  " + "  ".join(parts))

    allr = [v for d in per_day.values() for v in d["ratio"]]
    p = analysis.pool_ci(allr)
    if p.get("n"):
        print(f"\n全期斜率比 bare/enc: mean={p['mean']:.2f}× "
              f"CI[{p['lo']:.2f},{p['hi']:.2f}] median={p['median']:.2f} n={p['n']}")


if __name__ == "__main__":
    main()
