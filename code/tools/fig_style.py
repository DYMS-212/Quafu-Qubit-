"""统一图表样式 v2 —— Nature 系投稿级规范(全部报告图共用)。

图表契约(每图五点,定稿别随手改):
1. 每图先有一句话结论(写在各图工具的 METRICS 表 title 里);
2. 面板只留携带独立证据的;3. 原型=quantitative grid;
4. 后端=Python/matplotlib 独占;5. 导出=SVG(文字可编辑)+PDF(Type42)+PNG(300dpi)。

配色政策(克制,一图一套):
- 主角(编码/我们的方法)= 藏蓝 blue_main;第二主角(融合/相位)= teal;
- 基线(物理比特/对照)= 中性灰家族;红/绿只做方向性提示与死比特状态,绝不当系列色;
- 量级(热图)= 单色蓝渐变;文字永远用文字色。
- 单轴纪律;直接标注优先于图例;网格退后。

在 Mac 上跑(服务器无 CJK 字体)。
"""
from __future__ import annotations

import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

FIGROOT = Path(__file__).resolve().parent.parent / "figures"


def outdir() -> Path:
    """产物按生成日期归档:figures/<YYYY-MM-DD>/(同日重跑覆盖同名文件)。"""
    d = FIGROOT / time.strftime("%Y-%m-%d")
    d.mkdir(parents=True, exist_ok=True)
    return d

# ------------------------------------------------------------------ 角色化调色板
SURFACE = "#ffffff"          # 论文图一律纯白
INK = "#272727"              # 主文字(neutral_black)
INK2 = "#4D4D4D"             # 次文字
MUTED = "#767676"            # 弱文字/轴刻度
GRID = "#E8E7E3"             # 网格(极浅)
BASELINE = "#BDBBB6"         # 轴线

ROLE = {
    "hero":      "#0F4D92",  # 编码/我们的方法(blue_main)
    "hero2":     "#3775BA",  # 主角家族第二档(blue_secondary)
    "hero3":     "#86B6EF",  # 主角家族第三档(浅)
    "accent":    "#42949E",  # 第二主角:融合/相位(teal)
    "violet":    "#9A4D8E",  # 第三对照(解码器等多臂时)
    "base":      "#767676",  # 物理基线(neutral_mid)
    "base2":     "#CFCECE",  # 次级基线(neutral_light)
    "good":      "#2E9E44",  # 方向性:改善(专用,不当系列)
    "bad":       "#B64342",  # 方向性:恶化/死比特/失效(专用)
}
EMPH_GRAY = ROLE["base"]
CRITICAL = ROLE["bad"]

# 单色蓝渐变(量级;浅=小,深=大),与 hero 同家族
_SEQ_BLUE = ["#E8F1FC", "#CDE2FB", "#9EC5F4", "#6DA7EC", "#3987E5",
             "#2A78D6", "#1C5CAB", "#0F4D92", "#0D366B"]


def seq_cmap() -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list("seq_blue", _SEQ_BLUE)


def apply_style(font_size: float = 9):
    plt.rcParams.update({
        # 字体:Hiragino 在前保 CJK(本版 mpl 无逐字回退);
        # 对数轴一律用 log_axis() 纯文本刻度,避免 mathtext 负号缺字
        "font.sans-serif": ["Hiragino Sans GB", "PingFang SC",
                            "Arial Unicode MS", "Arial", "DejaVu Sans"],
        "font.family": "sans-serif",
        "axes.unicode_minus": False,
        "mathtext.fontset": "dejavusans",
        # 可编辑文字导出(nature-figure 强制规则)
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        # 版面
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "text.color": INK,
        "axes.edgecolor": BASELINE,
        "axes.labelcolor": INK2,
        "axes.titlecolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.axisbelow": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
        "lines.linewidth": 1.8,
        "font.size": font_size,
        "axes.titlesize": font_size + 1.5,
        "figure.dpi": 110,
    })


def panel_label(ax, label: str, x: float = -0.08, y: float = 1.04):
    """Nature 式面板字母(左上,粗体小写)。"""
    ax.text(x, y, label, transform=ax.transAxes, fontsize=12,
            fontweight="bold", color=INK, ha="left", va="bottom")


def log_axis(ax):
    """对数 y 轴 + 纯文本小数刻度(0.004/0.01/0.05…),不走 mathtext。"""
    from matplotlib.ticker import FuncFormatter, LogLocator
    ax.set_yscale("log")
    fmt = FuncFormatter(lambda v, _: f"{v:g}")
    ax.yaxis.set_major_formatter(fmt)
    ax.yaxis.set_minor_formatter(fmt)
    ax.yaxis.set_minor_locator(LogLocator(subs=(2, 3, 5)))


ERRKW = dict(capsize=2.5, elinewidth=1.1, capthick=1.1)


def save(fig, name: str, dpi: int = 300, formats=("png",)) -> list[Path]:
    """报告级图传 formats=("svg","pdf","png");讨论用默认只出 png。落到 outdir()。"""
    base = outdir() / name
    outs = []
    for fmt in formats:
        p = base.with_suffix(f".{fmt}")
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        outs.append(p)
    plt.close(fig)
    print("已写 " + ", ".join(str(p) for p in outs))
    return outs


def place_end_labels(ax, ends: list[tuple]):
    """线尾直接标注,像素空间防碰撞(自上而下强制最小间距)。
    ends: [(x, y, label, color)];须在坐标轴 scale/margins 定稿后调用。"""
    if not ends:
        return
    ax.autoscale_view()
    disp = sorted(((float(ax.transData.transform((x, y))[1]), x, y, lab, col)
                   for x, y, lab, col in ends), key=lambda t: -t[0])
    min_gap = 13.0
    placed = []
    for i, (py, x, y, lab, col) in enumerate(disp):
        if placed and placed[-1][0] - py < min_gap:
            py = placed[-1][0] - min_gap
        placed.append((py, x, y, lab, col))
    inv = ax.transData.inverted()
    for py, x, y, lab, col in placed:
        y_new = inv.transform((0, py))[1]
        ax.annotate(lab, (x, y_new), xytext=(7, 0),
                    textcoords="offset points", va="center",
                    fontsize=8, color=INK2)
