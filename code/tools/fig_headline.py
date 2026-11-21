"""headline 报告主图(零额度,吃 results/ 全历史,数据多一天图自动长一天)。

三张(报告级 svg+pdf+png,配同名 CSV 供核查;图上只说人话,代号只进文件名):
1. headline_e24_duel     —— 逻辑 vs 物理 Bell 保真度逐日对决 + 接受率随行面板
2. headline_e24_sentinel —— 接受率 × 逻辑保真度散点(接受率=免费的质量哨兵)
3. headline_e05_distance —— 重复码逻辑错误率 p_L(d) 距离扫描(逐日细线+跨日聚合)

统计口径(用户铁律):有效日自 07-08 起;跨日聚合 = 日均值序列的
mean ± bootstrap 95% CI(日间漂移是主不确定性);后选数字与接受率同图呈现。

运行(Mac): ./qrun tools/fig_headline.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))

import matplotlib.pyplot as plt  # noqa: E402

from qfab import analysis  # noqa: E402
from tools import fig_style as st  # noqa: E402

R = PROJ / "results"
FIRST_DAY = "2026-07-08"          # 07-07 首日队列大堵、数据薄,不入 headline 统计


def _date_of(stem: str) -> str:
    return next((p for p in stem.split("_") if len(p) == 10 and p[4] == "-"),
                stem[-10:])


def e24_rounds(chip: str = "Baihua") -> list[dict]:
    """每轮一条:date/logical/physical/acceptance(有效日起,按时间排序)。"""
    out = []
    for f in sorted((R / "e24_logical_bell").glob(f"summary_{chip}_*.json")):
        date = _date_of(f.stem)
        if date < FIRST_DAY:
            continue
        arms = json.loads(f.read_text()).get("arms", {})
        lo, ph = arms.get("logical") or {}, arms.get("physical") or {}
        if lo.get("fidelity") is None or ph.get("fidelity") is None:
            continue
        out.append(dict(date=date, run=f.stem, logical=lo["fidelity"],
                        physical=ph["fidelity"], acceptance=lo.get("acceptance")))
    return out


def _write_csv(path: Path, rows: list[dict]):
    if not rows:
        return
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"已写 {path}")


def _daily(rounds: list[dict], key: str) -> dict[str, dict]:
    by: dict[str, list[float]] = {}
    for r in rounds:
        if isinstance(r.get(key), (int, float)):
            by.setdefault(r["date"], []).append(r[key])
    return {d: analysis.pool_ci(v) for d, v in sorted(by.items())}


# ------------------------------------------------------------------ 图 1:对决
def fig_duel():
    rounds = e24_rounds()
    if not rounds:
        print("(无 e24 数据,跳过 duel)")
        return
    lo, ph, ac = (_daily(rounds, k) for k in ("logical", "physical", "acceptance"))
    days = sorted(lo)
    x = range(len(days))
    wins = sum(lo[d]["mean"] > ph[d]["mean"] for d in days)

    st.apply_style()
    fig, (a1, a2) = plt.subplots(
        2, 1, figsize=(7.0, 4.6), sharex=True,
        gridspec_kw=dict(height_ratios=[2.6, 1.0], hspace=0.12))

    for ax, series, color, lab in ((a1, lo, st.ROLE["hero"], "逻辑比特对"),
                                   (a1, ph, st.ROLE["base"], "物理比特对")):
        m = [series[d]["mean"] for d in days]
        err = ([series[d]["mean"] - series[d]["lo"] for d in days],
               [series[d]["hi"] - series[d]["mean"] for d in days])
        ax.errorbar(x, m, yerr=err, color=color, marker="o", ms=3.5,
                    **st.ERRKW, label=lab)
    for i, d in enumerate(days):                      # 失利日:方向性红点标注
        if lo[d]["mean"] <= ph[d]["mean"]:
            a1.plot(i, lo[d]["mean"], "o", ms=6, mfc="none",
                    mec=st.ROLE["bad"], mew=1.4)
    a1.set_ylabel("Bell 保真度(检错后选)")
    st.place_end_labels(a1, [
        (len(days) - 1, lo[days[-1]]["mean"], "逻辑比特对", st.ROLE["hero"]),
        (len(days) - 1, ph[days[-1]]["mean"], "物理比特对", st.ROLE["base"])])
    a1.set_title(f"编码层纠缠 {wins}/{len(days)} 个有效日高于物理层", loc="left")
    a1.set_title("红圈=失利日;点=日内聚合±95%CI", loc="right",
                 fontsize=8, color=st.MUTED)
    st.panel_label(a1, "a")

    am = [ac[d]["mean"] if d in ac else None for d in days]
    a2.plot(x, am, color=st.ROLE["accent"], marker="o", ms=3.5)
    a2.set_ylabel("接受率")
    a2.set_ylim(0.4, 0.85)
    for i, d in enumerate(days):
        if d in ac and lo[d]["mean"] <= ph[d]["mean"]:
            a2.plot(i, ac[d]["mean"], "o", ms=6, mfc="none",
                    mec=st.ROLE["bad"], mew=1.4)
    a2.set_title("失利日均伴随接受率骤降(同日红圈)", loc="left",
                 fontsize=9, color=st.INK2)
    a2.set_title("缺 14/15 日=平台故障", loc="right",
                 fontsize=8, color=st.MUTED)
    st.panel_label(a2, "b")
    # 全部数据都在 7 月,横轴只标日号(用户定)
    a2.set_xticks(list(x), [str(int(d[8:])) for d in days])
    a2.set_xlabel("7 月 · 日")

    st.save(fig, "headline_e24_duel", formats=("svg", "pdf", "png"))
    _write_csv(st.outdir() / "headline_e24_duel.csv", [
        dict(date=d, logical_mean=round(lo[d]["mean"], 5),
             logical_lo=round(lo[d]["lo"], 5), logical_hi=round(lo[d]["hi"], 5),
             physical_mean=round(ph[d]["mean"], 5),
             acceptance=round(ac[d]["mean"], 5) if d in ac else None,
             n_rounds=lo[d]["n"]) for d in days])


# ------------------------------------------------------------------ 图 2:哨兵
def fig_sentinel():
    rounds = [r for r in e24_rounds() if isinstance(r.get("acceptance"), float)]
    if len(rounds) < 3:
        print("(e24 带接受率的轮次不足,跳过 sentinel)")
        return
    win = [r for r in rounds if r["logical"] > r["physical"]]
    loss = [r for r in rounds if r["logical"] <= r["physical"]]
    xs = [r["acceptance"] for r in rounds]
    ys = [r["logical"] for r in rounds]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    r_pear = sxy / (sxx * syy) ** 0.5 if sxx and syy else float("nan")

    st.apply_style()
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    ax.scatter([r["acceptance"] for r in win], [r["logical"] for r in win],
               s=26, color=st.ROLE["hero"], label="当轮逻辑>物理")
    ax.scatter([r["acceptance"] for r in loss], [r["logical"] for r in loss],
               s=30, color=st.ROLE["bad"], marker="D", label="当轮逻辑≤物理")
    if sxx:
        k = sxy / sxx
        xf = [min(xs), max(xs)]
        ax.plot(xf, [my + k * (v - mx) for v in xf], ls="--", lw=1.1,
                color=st.ROLE["base2"], zorder=0)
    ax.set_xlabel("检错后选接受率")
    ax.set_ylabel("逻辑 Bell 保真度")
    ax.set_title("接受率下降,预警逻辑纠缠质量下降", loc="left")
    ax.set_title(f"每点=一轮;r={r_pear:.2f}", loc="right",
                 fontsize=8, color=st.MUTED)
    ax.legend(loc="lower right", fontsize=8)
    st.save(fig, "headline_e24_sentinel", formats=("svg", "pdf", "png"))
    _write_csv(st.outdir() / "headline_e24_sentinel.csv", [
        dict(date=r["date"], run=r["run"], acceptance=r["acceptance"],
             logical=r["logical"], physical=r["physical"],
             wins=r["logical"] > r["physical"]) for r in rounds])


# ------------------------------------------------------------------ 图 3:距离扫描
def fig_distance():
    per_day: dict[str, dict[int, list[float]]] = {}
    for f in sorted((R / "e05_repcode").glob("summary_2*.json")):
        date = _date_of(f.stem)
        if date < FIRST_DAY:
            continue
        data = json.loads(f.read_text())
        rows = data.get("rows", data) if isinstance(data, dict) else data
        for r in rows:
            if r.get("code") == "bit" and r.get("bit") == 1 and "err" in r:
                per_day.setdefault(date, {}).setdefault(r["d"], []).append(r["err"])
    if not per_day:
        print("(无 e05 数据,跳过 distance)")
        return
    ds = sorted({d for by in per_day.values() for d in by})
    daily_mean = {date: {d: sum(v) / len(v) for d, v in by.items()}
                  for date, by in per_day.items()}
    pooled = {d: analysis.pool_ci([m[d] for m in daily_mean.values() if d in m])
              for d in ds}

    st.apply_style()
    fig, ax = plt.subplots(figsize=(4.8, 3.6))
    for date, m in sorted(daily_mean.items()):        # 逐日细线(轻)
        dd = [d for d in ds if d in m]
        ax.plot(dd, [m[d] for d in dd], color=st.ROLE["hero3"],
                lw=0.8, alpha=0.55, zorder=1)
    mm = [pooled[d]["mean"] for d in ds]
    err = ([pooled[d]["mean"] - pooled[d]["lo"] for d in ds],
           [pooled[d]["hi"] - pooled[d]["mean"] for d in ds])
    ax.errorbar(ds, mm, yerr=err, color=st.ROLE["hero"], marker="o",
                ms=4.5, lw=2.0, **st.ERRKW, zorder=3)
    st.log_axis(ax)
    ax.set_xticks(ds)
    ax.set_xlabel("码距 d")
    ax.set_ylabel("逻辑错误率(多数表决)")
    n_days = len(daily_mean)
    # 诚实标题:数据显示 d≥9 后触地板(关联噪声/SPAM 下限),不吹"一直压"
    ax.set_title("加大码距持续压低错误率,至 d≈9 触地板", loc="left")
    ax.text(0.985, 0.97, f"粗=跨{n_days}日均值±95%CI\n细=单日",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=7.5, color=st.MUTED)
    st.save(fig, "headline_e05_distance", formats=("svg", "pdf", "png"))
    _write_csv(st.outdir() / "headline_e05_distance.csv", [
        dict(d=d, mean=round(pooled[d]["mean"], 6), lo=round(pooled[d]["lo"], 6),
             hi=round(pooled[d]["hi"], 6), n_days=pooled[d]["n"]) for d in ds])


if __name__ == "__main__":
    fig_duel()
    fig_sentinel()
    fig_distance()
