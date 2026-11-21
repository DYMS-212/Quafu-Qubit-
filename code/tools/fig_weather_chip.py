"""芯片"天气图":全芯片读出误差逐日/逐次体检的空间分布(零额度,吃 e14 数据)。

一图三产:
  figures/weather_chip_latest.png   最新一次体检的静态图(报告用)
  figures/weather_chip_days.png     小倍图:每天早体检一格(报告用,漂移一眼可见)
  figures/weather_chip.gif          时间流逝动画(答辩用,逐帧=每次体检)

编码:格色 = 读出误差 1-η(蓝色渐变,深=差;色标全帧固定,可跨帧比较);
死比特(当日不可用)= 红 ×。布局 = 芯片真实 12×13 网格(q = 行*13 + 列)。

在 Mac 上跑: ./qrun tools/fig_weather_chip.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import animation, colors  # noqa: E402

from qfab import analysis  # noqa: E402
from tools import fig_style as st  # noqa: E402

ROWS, COLS = 12, 13
VMAX = 0.5      # 1-η ≥ 0.5 已是废比特量级,截断保证色标分辨率
PROBE = PROJ / "results" / "e14_probe"


def frames() -> list[dict]:
    """全部体检帧:[{date, time, eta: {q: eta}}],按时间排序。"""
    out = []
    for day in sorted(d for d in PROBE.iterdir() if d.is_dir()):
        f0s = sorted(day.glob("selfcal_prep0__*.json"))
        f1s = sorted(day.glob("selfcal_prep1__*.json"))
        for f0, f1 in zip(f0s, f1s):
            r0, r1 = json.loads(f0.read_text()), json.loads(f1.read_text())
            if not (r0.get("result", {}).get("count")
                    and r1.get("result", {}).get("count")):
                continue
            qubits = r0["meta"]["qubits"]
            asg = analysis.assignment_from_p1(
                analysis.bit_p1(r0["result"]["count"], len(qubits)),
                analysis.bit_p1(r1["result"]["count"], len(qubits)))
            ts = (r0.get("submitted_at") or "")[11:16]
            out.append({"date": day.name, "time": ts,
                        "eta": {q: asg[i]["eta"] for i, q in enumerate(qubits)}})
    return out


def grid_of(eta: dict) -> np.ndarray:
    g = np.full((ROWS, COLS), np.nan)
    for q, e in eta.items():
        g[q // COLS, q % COLS] = 1 - e
    return g


def draw(ax, fr: dict, show_ids: bool = False):
    g = grid_of(fr["eta"])
    cmap = st.seq_cmap()
    cmap.set_bad(st.SURFACE)
    im = ax.pcolormesh(np.clip(g, 0, VMAX), cmap=cmap,
                       vmin=0, vmax=VMAX, edgecolors=st.SURFACE, linewidth=1.2)
    for q in range(ROWS * COLS):
        r, c = q // COLS, q % COLS
        if q not in fr["eta"]:
            ax.text(c + 0.5, r + 0.5, "×", ha="center", va="center",
                    color=st.CRITICAL, fontsize=9, fontweight="bold")
        elif show_ids:
            err = 1 - fr["eta"][q]
            ink = "#ffffff" if err > 0.55 * VMAX else st.MUTED
            ax.text(c + 0.5, r + 0.5, str(q), ha="center", va="center",
                    color=ink, fontsize=4.5)
    ax.set_xlim(0, COLS)
    ax.set_ylim(ROWS, 0)               # 行 0 在上
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    return im


def stats_line(fr: dict) -> str:
    errs = [1 - e for e in fr["eta"].values()]
    dead = ROWS * COLS - len(fr["eta"])
    return (f"误差中位 {np.median(errs):.3f} · 差比特(1-η>0.3) "
            f"{sum(e > 0.3 for e in errs)} 个 · 不可用 {dead} 个")


def main():
    st.apply_style()
    frs = frames()
    if not frs:
        sys.exit("无 e14 体检数据")
    print(f"共 {len(frs)} 帧体检数据({frs[0]['date']} → {frs[-1]['date']})")

    # 1. 最新静态图(带比特号)
    fr = frs[-1]
    fig, ax = plt.subplots(figsize=(7.4, 6.2))
    im = draw(ax, fr, show_ids=True)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label(f"读出误差 1-η(≥{VMAX} 截断)", color=st.INK2)
    cb.outline.set_visible(False)
    ax.set_title(f"芯片读出误差天气图 · Baihua · {fr['date']} {fr['time']}(自测)",
                 pad=10)
    ax.text(0, ROWS + 0.55, stats_line(fr) + " · × = 当日不可用",
            fontsize=8.5, color=st.MUTED)
    st.save(fig, "weather_chip_latest", formats=("png", "svg", "pdf"))

    # 2. 小倍图:每天第一次(早)体检
    days: dict[str, dict] = {}
    for f in frs:
        days.setdefault(f["date"], f)
    dl = list(days.values())
    ncol = min(4, len(dl))
    nrow = (len(dl) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 2.9 * nrow + 0.5),
                             squeeze=False)
    for ax in axes.flat:
        ax.set_visible(False)
    for ax, f in zip(axes.flat, dl):
        ax.set_visible(True)
        im = draw(ax, f)
        errs = [1 - e for e in f["eta"].values()]
        ax.set_title(f"{int(f['date'][8:])}日 · 中位 {np.median(errs):.3f}",
                     fontsize=9.5, pad=4)
    cb = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02)
    cb.set_label("读出误差 1-η", color=st.INK2)
    cb.outline.set_visible(False)
    fig.suptitle("同一块芯片,每天都是新天气 —— 每日自测读出误差", y=0.99)
    st.save(fig, "weather_chip_days", formats=("png", "svg", "pdf"))

    # 3. GIF 动画(逐帧=每次体检,含日内早晚)
    fig, ax = plt.subplots(figsize=(6.6, 5.9))
    im = draw(ax, frs[0])
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label("读出误差 1-η", color=st.INK2)
    cb.outline.set_visible(False)
    title = ax.set_title("", pad=10)
    sub = ax.text(0, ROWS + 0.55, "", fontsize=8.5, color=st.MUTED)

    def update(i):
        ax.collections[0].remove()
        for t in list(ax.texts):
            if t not in (sub,):
                t.remove()
        draw(ax, frs[i])
        title.set_text(f"芯片天气 · {frs[i]['date']} {frs[i]['time']}")
        sub.set_text(stats_line(frs[i]))
        return []

    ani = animation.FuncAnimation(fig, update, frames=len(frs))
    out = st.outdir() / "weather_chip.gif"
    ani.save(out, writer=animation.PillowWriter(fps=1), dpi=130)
    plt.close(fig)
    print(f"已写 {out}")

    # 4. 每帧统计 CSV(图配数,可核查;逐比特原始 η 见 data_export/<日期>/e14_chip_eta.csv)
    import csv
    p = st.outdir() / "weather_chip_stats.csv"
    with p.open("w", newline="") as fcsv:
        w = csv.DictWriter(fcsv, fieldnames=["date", "time", "err_median",
                                             "n_bad_gt0.3", "n_dead"])
        w.writeheader()
        for f in frs:
            errs = [1 - e for e in f["eta"].values()]
            w.writerow({"date": f["date"], "time": f["time"],
                        "err_median": round(float(np.median(errs)), 4),
                        "n_bad_gt0.3": sum(e > 0.3 for e in errs),
                        "n_dead": ROWS * COLS - len(f["eta"])})
    print(f"已写 {p}")


if __name__ == "__main__":
    main()
