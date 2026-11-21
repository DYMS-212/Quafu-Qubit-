"""全实验"气象统计":每个实验的关键量按天呈现(零额度,吃 results/ 全历史)。

每个实验一张单图(figures/weather_<exp>.png,讨论用)+ 四张支柱总览页
(figures/pillar_P*.png/svg/pdf,报告级,带面板字母)。
统一形式:当日逐轮散点(轻)+ 日均值±bootstrap 95% CI(重)+ 最佳/最差日标注;
配色=角色化(编码蓝主角/融合与相位 teal/物理灰基线);误差率对数轴;参考线虚灰。

运行(Mac): ./qrun tools/fig_weather_metrics.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))

import matplotlib.pyplot as plt  # noqa: E402

from qfab import analysis  # noqa: E402
from tools import fig_style as st  # noqa: E402

R = PROJ / "results"


def summaries(exp: str, pattern: str = "summary_*.json"):
    d = R / exp
    out = []
    if d.exists():
        for f in sorted(d.glob(pattern)):
            data = json.loads(f.read_text())
            rows = data.get("rows", data) if isinstance(data, dict) else data
            date = next((p for p in f.stem.split("_") if len(p) == 10
                         and p[4] == "-"), f.stem[-10:])
            out.append((date, rows, data))
    return out


def collect(exp: str, extractor, pattern: str = "summary_2*.json"):
    by_day: dict[str, list[float]] = {}
    for date, rows, data in summaries(exp, pattern):
        try:
            v = extractor(rows, data)
        except Exception:  # noqa: BLE001
            v = None
        if isinstance(v, (int, float)):
            by_day.setdefault(date, []).append(float(v))
    return by_day


# ------------------------------------------------------------------ 提取器
def first(it):
    return next(iter(it), None)


def x_e02_eta(rows, _):
    etas = sorted(r["eta"] for r in rows)
    return etas[len(etas) // 2]


def x_e03(field):
    def f(rows, _):
        vals = sorted(abs(r["asym"]) if field == "asym" else r["bell_mean"]
                      for r in rows)
        return vals[len(vals) // 2]
    return f


def x_e17(kind):
    return lambda rows, data: data.get(kind, {}).get("excess_ratio")


def x_e08(rows, _):
    vals = sorted(r["delta"] for r in rows if not r.get("inject"))
    return vals[len(vals) // 2] if vals else None


def x_e05(d):
    return lambda rows, _: first(r["err"] for r in rows
                                 if r.get("code") == "bit" and r.get("d") == d
                                 and r.get("bit") == 1)


def x_e10_star(rows, _):
    stars: dict = {}
    for r in rows:
        if "ps_err_perq" in r and "kind" not in r:
            stars.setdefault(r["star"], []).append(r["ps_err_perq"])
    return min((sum(v) / len(v) for v in stars.values()), default=None)


def x_e10_fusion(rows, _):
    paths: dict = {}
    for r in rows:
        if r.get("kind") == "fusion":
            paths.setdefault(r["vq"], []).append(r["ps_err_perq"])
    return min((sum(v) / len(v) for v in paths.values()), default=None)


def x_e10_phys(rows, _):
    per_q: dict = {}
    for r in rows:
        if r.get("kind") == "phys":
            per_q.setdefault(r["qubit"], []).append(r["err"])
    return min((sum(v) / len(v) for v in per_q.values()), default=None)


def x_e11(code, kind):
    def f(rows, _):
        vals = [r["err"] for r in rows if r.get("kind") == kind
                and r.get("bit") == 1 and r.get("K") == 32
                and r.get("code", "bit") == code
                and (kind == "phys" or r.get("d") == 7)]
        return min(vals) if vals else None
    return f


def x_e19_slope(field):
    return lambda rows, _: analysis.fit_decay_per_layer(
        [(r["k"], r.get(field)) for r in rows]).get("r")


def x_e20(rows, _):
    fro = first(r["err"] for r in rows if r.get("kind") == "rep"
                and r.get("d") == 7 and r.get("arm") == "frozen")
    day = first(r["err"] for r in rows if r.get("kind") == "rep"
                and r.get("d") == 7 and r.get("arm") == "daily")
    return fro / day if fro and day else None


def x_e21(metric):
    # 两指标选中同一条链时实验会合并成一行 metrics=['readout','composite'],用包含匹配
    return lambda rows, _: first(r["err"] for r in rows
                                 if metric in (r.get("metrics") or []))


def x_e06(field):
    return lambda rows, _: first(r[field] for r in rows
                                 if r.get("code") == "bit" and r.get("d") == 7
                                 and r.get("bit") == 1)


def x_e12(field):
    return lambda rows, _: first(r.get(field) for r in rows
                                 if r.get("code") == "bit" and r.get("d") == 7
                                 and r.get("bit") == 1)


def x_e13(kind):
    def f(rows, _):
        by = {}
        for r in rows:
            by.setdefault(tuple(r.get("edge", [])), {})[r.get("basis")] = r
        for pair in by.values():
            if "ZZ" in pair and "XX" in pair:
                return abs(pair["ZZ"].get(kind, 0) - pair["XX"].get(kind, 0))
        return None
    return f


def x_e04(rows, _):
    vals = [r["success"] for r in rows if "success" in r]
    return sum(vals) / len(vals) if vals else None


def x_e07(rows, _):
    p = {r["case"]: r["p_flag_1"] for r in rows}
    return p.get("X_anc", 0) - p.get("none", 0) if "X_anc" in p else None


def x_e09(layout, n=13):
    """钉死 N=13:跨日可比;当日没跑到 N=13 就空着(诚实),不偷换成小 N。"""
    return lambda rows, _: first(
        r["value"] for r in rows if r.get("layout") == layout
        and r.get("N") == n and r.get("metric") == "pop00_11")


def x_e18(which):
    def f(rows, _):
        if which == "none":
            return first(r["p_expected"] for r in rows if r["case"] == "none")
        if which == "mid":
            vals = [r["p_expected"] for r in rows if "_mid" in r["case"]]
            return min(vals) if vals else None
        vals = [r["p_expected"] for r in rows if "pre" in r["case"]]
        return first(vals)
    return f


def x_e22(arm):
    def f(rows, _):
        vals = [r["err"] for r in rows if r.get("arm") == arm]
        return sum(vals) / len(vals) if vals else None
    return f


def x_e25(kind):
    """码族规模轴:k=8, K=8, L1 的每逻辑比特错误率(enc vs bare)。"""
    def f(rows, _):
        vals = [r["err_per_logical"] for r in rows
                if r.get("kind") == kind and r.get("k") == 8 and r.get("K") == 8
                and r.get("logical") == 1 and r.get("err_per_logical") is not None]
        return vals[0] if vals else None
    return f


def x_e26(with_t):
    def f(rows, _):
        vals = [r["corr"] for r in rows if r.get("kind") == "ccz"
                and r.get("with_t") == with_t and r.get("corr") is not None]
        return sum(vals) / len(vals) if vals else None
    return f


def x_e24(arm_):
    return lambda rows, data: (data.get("arms") or {}).get(arm_, {}).get("fidelity")


def x_e23_bell(rows, data):
    return (data.get("states") or {}).get("bell", {}).get("fidelity")


# ------------------------------------------------------------------ 图表定义
def arm(label, exp, ex, color, pattern="summary_2*.json", agg="mean"):
    return dict(label=label, exp=exp, ex=ex, color=color, pattern=pattern, agg=agg)


H, H2, H3, A, V, B, B2 = (st.ROLE["hero"], st.ROLE["hero2"], st.ROLE["hero3"],
                          st.ROLE["accent"], st.ROLE["violet"],
                          st.ROLE["base"], st.ROLE["base2"])

# 图上一律说人话:实验代号只出现在文件名与代码里,不出现在标题/标注
PILLARS = [
    ("P1_device", "支柱一 · 设备真相与统计自证", [
        dict(better="high", name="e02_eta", title="全芯片读出质量 η(每日自测,中位)", log=False,
             arms=[arm("η 中位", "e02_selfcal", x_e02_eta, H, "assignment_*.json")]),
        dict(better="high", name="e03_bell", title="双比特纠缠边质量(全芯片中位)", log=False,
             arms=[arm("Bell 保真代理", "e03_bell", x_e03("bell"), H, "atlas_Baihua_*.json"),
                   arm("相干不对称 |ZZ-XX|", "e03_bell", x_e03("asym"), A,
                       "atlas_Baihua_*.json")]),
        dict(name="e17_var", title="误差棒可信度:批内涨落 ÷ 散粒噪声(1 = 理想)",
             log=False, ref=1.0,
             arms=[arm("纠缠对", "e17_variance", x_e17("bell"), H),
                   arm("重复码", "e17_variance", x_e17("rep"), A)]),
        dict(name="e08_tiling", title="70 比特并行运行的串扰(中位增量)",
             log=False, ref=0.0,
             arms=[arm("并行 - 单独", "e08_tiling", x_e08, H)]),
    ]),
    ("P2_codes", "支柱二 · 编码对决(存储 / 静态 / 门级 / 选址)", [
        dict(name="e05_dist", title="码距越大压得越狠:重复码逻辑错误率", log=True,
             arms=[arm("码距 7", "e05_repcode", x_e05(7), H),
                   arm("码距 5", "e05_repcode", x_e05(5), H2),
                   arm("码距 3", "e05_repcode", x_e05(3), H3)]),
        dict(name="e10_duel", title="逻辑比特 vs 最强物理比特(静态存读)", log=True,
             arms=[arm("编码·星形制备", "e10_code422", x_e10_star, H),
                   arm("编码·验证制备", "e10_code422", x_e10_fusion, A),
                   arm("最强物理比特", "e10_code422", x_e10_phys, B)]),
        dict(name="e11_memory", title="编码存储:存满负载后的错误率(码距 7)", log=True,
             arms=[arm("编码·幅值(存|1>)", "e11_memory", x_e11("bit", "enc"), H),
                   arm("编码·相位(存|->)", "e11_memory", x_e11("phase", "enc"), A),
                   arm("物理·幅值", "e11_memory", x_e11("bit", "phys"), B),
                   arm("物理·相位", "e11_memory", x_e11("phase", "phys"), B2)]),
        dict(name="e19_slope", title="逻辑门对决:每层门的错误率(斜率)", log=False,
             arms=[arm("编码比特", "e19_gates", x_e19_slope("ps_err_perq"), H),
                   arm("最强物理比特", "e19_gates", x_e19_slope("bare_err_perq"), B)]),
        dict(better="high", name="e20_gain", title="当日自适应选址的增益(同批对照,倍数)",
             log=False, ref=1.0,
             arms=[arm("冻结链 ÷ 当日链", "e20_anchor_ab", x_e20, H)]),
        dict(better="high", name="e24_lbell",
             title="逻辑纠缠 vs 物理纠缠:Bell 保真度(检错后选)", log=False,
             arms=[arm("逻辑(Baihua)", "e24_logical_bell", x_e24("logical"), H,
                       "summary_Baihua_*.json"),
                   arm("物理(Baihua)", "e24_logical_bell", x_e24("physical"), B,
                       "summary_Baihua_*.json"),
                   arm("逻辑(Shenglian)", "e24_logical_bell", x_e24("logical"), A,
                       "summary_Shenglian_*.json"),
                   arm("逻辑(Dongling)", "e24_logical_bell", x_e24("logical"), B2,
                       "summary_Dongling_*.json")]),
        dict(name="e25_iceberg", title="检错码族做大块:每逻辑比特错误率(k=8,存 8 层)", log=True,
             arms=[arm("编码 [[10,8,2]]", "e25_iceberg", x_e25("ice"), H),
                   arm("裸比特 ×8", "e25_iceberg", x_e25("bare"), B)]),
        dict(better="high", name="e26_ccz",
             title="逻辑 CCZ(magic 门):hypergraph 相关子 vs 无 T 对照", log=False,
             arms=[arm("打 T 层(理想 1.0)", "e26_code832", x_e26(True), H),
                   arm("无 T 对照(理想 0.5)", "e26_code832", x_e26(False), B)]),
        dict(name="e21_sel", title="选址靠什么:自测读出 vs 复合指标", log=True,
             arms=[arm("自测读出 η", "e21_selection", x_e21("readout"), H),
                   arm("η×厂商CZ 复合", "e21_selection", x_e21("composite"), V)]),
    ]),
    ("P3_decode", "支柱三 · 解码 · 后选 · 缓解", [
        dict(name="e06_ps", title="检错后选:错误率下降与接受率代价(码距 7)", log=True,
             arms=[arm("后选", "e06_postselect", x_e06("ps_err"), H),
                   arm("仅表决", "e06_postselect", x_e06("raw_majority_err"), B),
                   arm("接受率", "e06_postselect", x_e06("acceptance"), B2)]),
        dict(name="e12_decoder", title="解码器对决:朴素 MWPM 有害,最优≈表决(码距 7)",
             log=True,
             arms=[arm("精确 MAP", "e12_decoder", x_e12("G_map_qfit"), H),
                   arm("多数表决", "e12_decoder", x_e12("A_majority"), B),
                   arm("朴素 MWPM", "e12_decoder", x_e12("C_mwpm_uniform"), V)]),
        dict(name="e13_twirl", title="随机编译把相干误差搅成随机(不对称→0)",
             log=False, ref=0.0,
             arms=[arm("原始 |ZZ-XX|", "e13_mitigation", x_e13("raw"), B),
                   arm("随机编译后", "e13_mitigation", x_e13("twirl"), H)]),
        dict(better="high", name="e07_flag", title="flag 比特:传播错误的检出对比度", log=False,
             arms=[arm("检出 - 误报", "e07_flag", x_e07, H)]),
        dict(better="high", name="e04_synd", title="稳定子检验注错查找表:平均命中率", log=False,
             arms=[arm("命中率", "e04_syndrome", x_e04, H)]),
    ]),
    ("P4_scale", "支柱四 · 规模 · 时间维度 · 跨机", [
        dict(better="high", name="e09_ghz", title="会读拓扑的差距:13 比特 GHZ 布居", log=True,
             arms=[arm("树形布线", "e09_ghz", x_e09("tree"), H),
                   arm("线性布线", "e09_ghz", x_e09("linear"), B)]),
        dict(better="high", name="e18_t2", title="两轮稳定子:错误发生时刻的定位特征", log=False,
             arms=[arm("轮间错特征", "e18_t2swap", x_e18("mid"), H),
                   arm("无错基线", "e18_t2swap", x_e18("none"), B),
                   arm("编码前错特征", "e18_t2swap", x_e18("pre"), A)]),
        dict(better="high", name="e23_tomo",
             title="Bell 态层析保真度:三台机器(当日)", log=False,
             arms=[arm("Baihua", "e23_tomography", x_e23_bell, H,
                       "summary_Baihua_*.json"),
                   arm("Shenglian", "e23_tomography", x_e23_bell, A,
                       "summary_Shenglian_*.json"),
                   arm("Dongling", "e23_tomography", x_e23_bell, B,
                       "summary_Dongling_*.json")]),
        dict(name="e22_sigma", title="编码不压波动:轮间 σ 对比(诚实负结果)",
             log=False,
             arms=[arm("编码比特 σ", "e22_variance_panel", x_e22("enc"), H, "summary_*.json", "std"),
                   arm("裸比特 σ", "e22_variance_panel", x_e22("bare"), B, "summary_*.json", "std")]),
    ]),
]


# ------------------------------------------------------------------ 绘制
def plot_metric(ax, m, csv_rows: list | None = None):
    """零遮挡纪律:最佳/最差日放标题行右侧(图外),线尾标注防碰撞排布。
    csv_rows 非 None 时,把画的数字逐行填进去(每图配同名 CSV,可核查)。"""
    all_dates = sorted({d for a in m["arms"]
                        for d in collect(a["exp"], a["ex"], a["pattern"])})
    if not all_dates:
        return False
    xi = {d: i for i, d in enumerate(all_dates)}
    ends, hero = [], None
    for a in m["arms"]:
        by_day = collect(a["exp"], a["ex"], a["pattern"])
        xs, means, los, his, dates = [], [], [], [], []
        for d in all_dates:
            vals = by_day.get(d)
            if not vals:
                continue
            p = analysis.pool_ci(vals)
            if a["agg"] == "std":
                if p.get("n", 0) < 2:
                    continue
                v = p["sd"]
                xs.append(xi[d]); means.append(v)
                los.append(0); his.append(0); dates.append(d)
            else:
                xs.append(xi[d]); means.append(p["mean"])
                los.append(p["mean"] - p["lo"]); his.append(p["hi"] - p["mean"])
                dates.append(d)
                jit = [xi[d] + (j - len(vals) / 2) * 0.5 / max(len(vals), 1)
                       for j in range(len(vals))]
                ax.scatter(jit, vals, s=10, color=a["color"], alpha=0.3,
                           linewidths=0, zorder=2)
            if csv_rows is not None:
                csv_rows.append({"arm": a["label"], "date": d,
                                 "agg": a["agg"], "value": means[-1],
                                 "ci_lo": round(means[-1] - los[-1], 6),
                                 "ci_hi": round(means[-1] + his[-1], 6),
                                 "n_rounds": len(vals),
                                 "rounds": json.dumps([round(v, 6) for v in vals])})
        if not xs:
            continue
        ax.errorbar(xs, means, yerr=[los, his], color=a["color"], marker="o",
                    markersize=4.5, linewidth=1.8, zorder=3, **st.ERRKW)
        ends.append((xs[-1], means[-1], a["label"], a["color"]))
        if hero is None:
            hero = (dates, means)
    if m.get("ref") is not None:
        ax.axhline(m["ref"], color=st.BASELINE, linestyle="--",
                   linewidth=1.0, zorder=1)
    if m["log"]:
        st.log_axis(ax)
    ax.set_xticks(range(len(all_dates)))
    # 全部数据都在 7 月,横轴只标日号(用户定;月份进页标题)
    ax.set_xticklabels([str(int(d[8:])) for d in all_dates])
    ax.set_title(m["title"], loc="left")
    ax.margins(x=0.18)
    st.place_end_labels(ax, ends)
    if hero and len(hero[1]) > 1:                      # 图外标注,永不遮挡数据
        dates, means = hero
        if m.get("better") == "high":                  # η/保真/命中率类:越大越好
            bi, wi = means.index(max(means)), means.index(min(means))
        else:                                          # 错误率类:越小越好
            bi, wi = means.index(min(means)), means.index(max(means))
        ax.text(1.0, 1.02, f"最佳 {means[bi]:.3g}({int(dates[bi][8:])}日) · "
                           f"最差 {means[wi]:.3g}({int(dates[wi][8:])}日)",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=7.5, color=st.MUTED)
    return True


def write_metric_csv(name: str, rows: list):
    import csv as _csv
    p = st.outdir() / f"weather_{name}.csv"
    with p.open("w", newline="") as f:
        w = _csv.DictWriter(f, fieldnames=["arm", "date", "agg", "value",
                                           "ci_lo", "ci_hi", "n_rounds", "rounds"])
        w.writeheader()
        w.writerows(rows)


def main():
    st.apply_style()
    n_single = 0
    for pid, ptitle, metrics in PILLARS:
        drawn = []
        for m in metrics:
            fig, ax = plt.subplots(figsize=(6.4, 3.4))
            rows: list = []
            if plot_metric(ax, m, csv_rows=rows):
                ax.set_xlabel("日期 · 淡点=单轮 · 实心=日均值±95%CI")
                st.save(fig, f"weather_{m['name']}")
                write_metric_csv(m["name"], rows)
                drawn.append(m)
                n_single += 1
            else:
                plt.close(fig)
                print(f"跳过 {m['name']}(无数据)")
        if not drawn:
            continue
        ncol = 2
        nrow = (len(drawn) + ncol - 1) // ncol
        fig, axes = plt.subplots(nrow, ncol,
                                 figsize=(7.2 * ncol, 3.4 * nrow), squeeze=False)
        for ax in axes.flat:
            ax.set_visible(False)
        for i, (ax, m) in enumerate(zip(axes.flat, drawn)):
            ax.set_visible(True)
            plot_metric(ax, m)
            st.panel_label(ax, chr(ord("a") + i))
        fig.suptitle(f"{ptitle} · 7 月逐日分布(淡点=单轮,实心=日均值±95%CI)",
                     y=1.0, fontsize=13)
        fig.tight_layout(pad=1.8)
        st.save(fig, f"pillar_{pid}", formats=("png", "svg", "pdf"))
    print(f"\n单图 {n_single} 张(各配同名 CSV)+ 支柱页 {len(PILLARS)} 张 → {st.outdir()}")


if __name__ == "__main__":
    main()
