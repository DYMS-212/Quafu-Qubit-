"""终局统计引擎(零额度):从 results/ 全量重算报告要用的每一个数字。

定位(附录 B 的执行器):正文/答辩引用的数字一律出自本引擎,不手抄——
冻结(07-31)后重跑一次,输出即终版;之前任何一天跑,输出即"截至该日"版。

输出:data_export/freeze_stats/stats.json(机器可读,全部数字+溯源)
     data_export/freeze_stats/stats.md  (人读数字表,直接贴报告用)

统计口径(用户铁律,附录 B 1):
- headline = 跨日均值 ± 跨日 bootstrap 95% CI,n 的正主是天数;
  单日多轮的紧 CI 只说明"日内测量精度",绝不冒充总不确定性;
- 后选数字必附接受率;e24 按接受率双态分层(附录 B 3d),绝不报裸均值;
- 有效日自 07-08 起(07-07 首日队列大堵、数据薄,不入统计)。

运行(Mac): ./qrun tools/freeze_stats.py [--cutoff 2026-07-31]
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))

from qfab import analysis  # noqa: E402

R = PROJ / "results"
OUT = PROJ / "data_export" / "freeze_stats"
FIRST_DAY = "2026-07-08"
CHIPS = ("Baihua", "Dongling", "Shenglian")

WARN: list[str] = []


def _date_of(stem: str) -> str:
    return next((p for p in stem.split("_") if len(p) == 10 and p[4] == "-"),
                stem[-10:])


def _load(exp: str, cutoff: str, chip: str | None = None) -> list[tuple[str, dict]]:
    """(date, doc) 列表,按日期升序,限 [FIRST_DAY, cutoff]。chip=None 收全部。"""
    patt = f"summary_{chip}_*.json" if chip else "summary*.json"
    out = []
    for f in sorted((R / exp).glob(patt)):
        d = _date_of(f.stem)
        if FIRST_DAY <= d <= cutoff:
            out.append((d, json.loads(f.read_text())))
    if not out:
        WARN.append(f"{exp}({chip or '全部'}): 区间内无 summary")
    return out


def _pool(vals):
    p = analysis.pool_ci([v for v in vals if isinstance(v, (int, float))])
    return {k: round(v, 6) if isinstance(v, float) else v for k, v in p.items()}


def _daily_then_pool(day_vals: dict[str, list[float]]) -> dict:
    """日内先均值,再跨日 pool(n=天数)——headline 统一口径。"""
    means = {d: sum(v) / len(v) for d, v in sorted(day_vals.items()) if v}
    out = _pool(list(means.values()))
    out["daily_means"] = {d: round(m, 5) for d, m in means.items()}
    return out


def _pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    return round(sxy / (sxx * syy) ** 0.5, 4) if sxx and syy else None


# ---------------------------------------------------------------- 工作量
def workload(cutoff: str) -> dict:
    days = {}
    for f in sorted((R / "quota_log").glob("2026-*.txt")):
        d = f.stem
        if FIRST_DAY <= d <= cutoff:
            days[d] = sum(1 for _ in f.open())
    n_summaries = sum(1 for _ in R.glob("e*/summary*.json"))
    return {
        "有效天数": len(days), "起止": [min(days), max(days)] if days else None,
        "总任务数": sum(days.values()), "日均任务": round(sum(days.values()) / max(len(days), 1)),
        "逐日任务数": days, "summary文件总数": n_summaries,
        "机器": list(CHIPS),
    }


# ---------------------------------------------------------------- e24
def e24(cutoff: str) -> dict:
    out = {}
    for chip in CHIPS:
        rounds = []
        for d, doc in _load("e24_logical_bell", cutoff, chip):
            arms = doc.get("arms", {})
            lo, ph = arms.get("logical") or {}, arms.get("physical") or {}
            if lo.get("fidelity") is None or ph.get("fidelity") is None:
                continue
            rounds.append(dict(date=d, L=lo["fidelity"], P=ph["fidelity"],
                               acc=lo.get("acceptance")))
        if not rounds:
            continue
        day_L = defaultdict(list); day_P = defaultdict(list); day_acc = defaultdict(list)
        for r in rounds:
            day_L[r["date"]].append(r["L"]); day_P[r["date"]].append(r["P"])
            if r["acc"] is not None:
                day_acc[r["date"]].append(r["acc"])
        Ld, Pd = _daily_then_pool(day_L), _daily_then_pool(day_P)
        wins = sum(Ld["daily_means"][d] > Pd["daily_means"][d] for d in Ld["daily_means"])
        res = {
            "轮数": len(rounds), "天数": len(day_L),
            "跨日_逻辑F": Ld, "跨日_物理F": Pd, "跨日_接受率": _daily_then_pool(day_acc),
            "日胜负": f"{wins}/{len(Ld['daily_means'])}",
        }
        accs = [r["acc"] for r in rounds if r["acc"] is not None]
        Ls = [r["L"] for r in rounds if r["acc"] is not None]
        res["哨兵_全轮r"] = _pearson(accs, Ls)
        dm = Ld["daily_means"]
        am = {d: sum(v) / len(v) for d, v in day_acc.items()}
        common = sorted(set(dm) & set(am))
        res["哨兵_日级r"] = _pearson([am[d] for d in common], [dm[d] for d in common])
        # 双态分层(附录 B 3d):接受率 1 维 2-means,阈值=两簇中心中点
        if len(accs) >= 6:
            c1, c2 = min(accs), max(accs)
            for _ in range(50):
                g1 = [a for a in accs if abs(a - c1) <= abs(a - c2)]
                g2 = [a for a in accs if abs(a - c1) > abs(a - c2)]
                nc1 = sum(g1) / len(g1) if g1 else c1
                nc2 = sum(g2) / len(g2) if g2 else c2
                if abs(nc1 - c1) + abs(nc2 - c2) < 1e-9:
                    break
                c1, c2 = nc1, nc2
            thr = (c1 + c2) / 2
            good = [(r["L"], r["P"]) for r in rounds if (r["acc"] or 0) >= thr]
            bad = [(r["L"], r["P"]) for r in rounds if r["acc"] is not None and r["acc"] < thr]
            res["双态分层"] = {
                "阈值(2means)": round(thr, 4), "簇中心": [round(c1, 4), round(c2, 4)],
                "好态_逻辑F": _pool([g[0] for g in good]),
                "好态_物理F": _pool([g[1] for g in good]),
                "坏态_逻辑F": _pool([b[0] for b in bad]),
                "坏态_物理F": _pool([b[1] for b in bad]),
                "坏态占比": round(len(bad) / len(accs), 4),
                "好态轮数": len(good), "坏态轮数": len(bad),
            }
        out[chip] = res
    return out


# ---------------------------------------------------------------- e05
def e05(cutoff: str) -> dict:
    per = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))  # chip→code→d→[日内]
    for date, doc in _load("e05_repcode", cutoff):
        rows = doc.get("rows", doc) if isinstance(doc, dict) else doc
        for r in rows:
            if r.get("bit") == 1 and "err" in r:
                per[r.get("chip", "Baihua")][r["code"]][r["d"]].append((date, r["err"]))
    out = {}
    for chip, by_code in per.items():
        out[chip] = {}
        for code, by_d in by_code.items():
            curve = {}
            for d, pairs in sorted(by_d.items()):
                day_vals = defaultdict(list)
                for date, v in pairs:
                    day_vals[date].append(v)
                curve[f"d{d}"] = _daily_then_pool(day_vals)
            out[chip][code] = curve
        bit = out[chip].get("bit", {})
        if "d7" in bit and "d9" in bit:
            out[chip]["地板判定"] = {
                "d7_mean": bit["d7"]["mean"], "d9_mean": bit["d9"]["mean"],
                "d9/d7": round(bit["d9"]["mean"] / bit["d7"]["mean"], 3) if bit["d7"]["mean"] else None,
                "说明": "d9/d7≈1 即触地板(关联噪声/SPAM 下限),口径=压制至 d≈9 触地板",
            }
    return out


# ---------------------------------------------------------------- e11
def e11(cutoff: str) -> dict:
    per = defaultdict(lambda: defaultdict(list))
    phys_day = defaultdict(lambda: defaultdict(list))  # (chip,code,K)→date→[各臂err]
    site_moves = {}
    for date, doc in _load("e11_memory", cutoff):
        prov = doc.get("provenance", {})
        chip = prov.get("chip", "Baihua")
        qkey = str(prov.get("qubits", ""))[:60]
        site_moves.setdefault(chip, []).append((date, qkey))
        for r in doc.get("rows", []):
            # 口径=衰减臂 bit=1(|1⟩/叠加态才随 K 退化;bit=0 是平凡臂),与 digest/e05 一致
            if "err" not in r or "K" not in r or r.get("bit") != 1:
                continue
            code = r.get("code", "bit")
            if r.get("kind") == "enc" and "d" in r:
                per[(chip, code, r["d"], r["K"])][date].append(r["err"])
            elif r.get("kind") == "phys":
                phys_day[(chip, code, r["K"])][date].append(r["err"])
    out = defaultdict(dict)
    for (chip, code, d, K), day_vals in sorted(per.items()):
        out[chip][f"{code}_enc_d{d}_K{K}"] = _daily_then_pool(day_vals)
    for (chip, code, K), day_vals in sorted(phys_day.items()):
        # 3g 口径:物理基线取当日最好臂(min,对我们更严);中位臂另存
        out[chip][f"{code}_physbest_K{K}"] = _daily_then_pool(
            {d: [min(v)] for d, v in day_vals.items() if v})
        out[chip][f"{code}_physmed_K{K}"] = _daily_then_pool(
            {d: [sorted(v)[len(v) // 2]] for d, v in day_vals.items() if v})
    # headline 对决:d7 K 最大档 enc vs 最好物理臂(bit 与 phase 各一)
    for chip in list(out):
        for code in ("bit", "phase"):
            keys = [k for k in out[chip] if k.startswith(f"{code}_enc_d7")]
            if not keys:
                continue
            kmax = max(keys, key=lambda k: int(k.split("K")[-1]))
            Ktag = kmax.split("_")[-1]
            pk = f"{code}_physbest_{Ktag}"
            if pk in out[chip]:
                out[chip][f"headline_d7对决_{code}"] = {
                    "K档": Ktag, "enc": out[chip][kmax], "phys_best": out[chip][pk],
                    "CI分离": out[chip][kmax]["hi"] < out[chip][pk]["lo"],
                }
    # 换址记录(3g 口径)
    moves = {}
    for chip, seq in site_moves.items():
        change = [seq[0]]
        for d, q in seq[1:]:
            if q != change[-1][1]:
                change.append((d, q))
        moves[chip] = [d for d, _ in change]
    out["换址日期"] = moves
    return dict(out)


# ---------------------------------------------------------------- e10
def e10(cutoff: str) -> dict:
    day_best_star = defaultdict(list)
    day_phys = defaultdict(list)
    day_fusion = defaultdict(list)
    day_fusion_acc = defaultdict(list)
    day_star_acc = defaultdict(list)
    for date, doc in _load("e10_code422", cutoff):
        rows = doc.get("rows", [])
        stars = defaultdict(list)
        phys_by_q = defaultdict(list)   # 与 enc 同口径:每 qubit 对基底×bit 全平均
        for r in rows:
            lab, kind = str(r.get("label", "")), str(r.get("kind", ""))
            if kind == "phys_ss":
                phys_by_q[r.get("qubit")].append(r["err"])
            elif "fus" in lab or "fus" in kind:
                if r.get("ps_err_perq") is not None:
                    day_fusion[date].append(r["ps_err_perq"])
                if r.get("acceptance") is not None:
                    day_fusion_acc[date].append(r["acceptance"])
            elif r.get("ps_err_perq") is not None and r.get("star") is not None:
                stars[r["star"]].append(r)
        if phys_by_q:
            day_phys[date].append(min(sum(v) / len(v) for v in phys_by_q.values()))
        if stars:
            per_star = {s: sum(x["ps_err_perq"] for x in rs) / len(rs)
                        for s, rs in stars.items()}
            best = min(per_star, key=per_star.get)
            day_best_star[date].append(per_star[best])
            day_star_acc[date].extend(x["acceptance"] for x in stars[best]
                                      if x.get("acceptance") is not None)
    out = {
        "最好星_后选perq": _daily_then_pool(day_best_star),
        "最好星_接受率": _daily_then_pool(day_star_acc),
        "最强物理基线_perq": _daily_then_pool(day_phys),
        "融合制备_perq": _daily_then_pool(day_fusion),
        "融合_接受率": _daily_then_pool(day_fusion_acc),
    }
    bs, ph = out["最好星_后选perq"], out["最强物理基线_perq"]
    if bs.get("n") and ph.get("n"):
        out["对决判定"] = {"编码胜": bs["hi"] < ph["lo"], "说明": "跨日 CI 分离才算稳定差距"}
    fu = out["融合制备_perq"]
    if fu.get("n") and bs.get("n"):
        out["融合判定"] = {"融合优于星形": fu["hi"] < bs["lo"]}
    return out


# ---------------------------------------------------------------- e19
def e19(cutoff: str) -> dict:
    day_enc, day_bare, day_ratio, cross_days = (defaultdict(list), defaultdict(list),
                                                defaultdict(list), [])
    for date, doc in _load("e19_gates", cutoff):
        sl = doc.get("slopes", {})
        if (sl.get("enc") or {}).get("r") is not None and (sl.get("bare") or {}).get("r") is not None:
            day_enc[date].append(sl["enc"]["r"])
            day_bare[date].append(sl["bare"]["r"])
            if sl.get("ratio_bare_over_enc") is not None:
                day_ratio[date].append(sl["ratio_bare_over_enc"])
        ks = [r["k"] for r in doc.get("rows", [])
              if r.get("ps_err_perq") is not None and r.get("bare_err_perq") is not None
              and r["ps_err_perq"] < r["bare_err_perq"]]
        if ks:
            cross_days.append((date, min(ks)))
    enc_p, bare_p = _daily_then_pool(day_enc), _daily_then_pool(day_bare)
    em, bm = enc_p.get("daily_means", {}), bare_p.get("daily_means", {})
    common = sorted(set(em) & set(bm))
    ratio_daymean = {d: [bm[d] / em[d]] for d in common if em[d] > 1e-6}
    return {
        "enc斜率_跨日": enc_p,
        "bare斜率_跨日": bare_p,
        "enc斜率更低的天数": f"{sum(em[d] < bm[d] for d in common)}/{len(common)}",
        "斜率比_跨日": _daily_then_pool(ratio_daymean),
        "斜率比_注": "口径=日均bare斜率/日均enc斜率(轮级先除会被 enc≈0 轮爆炸污染,弃用)",
        "条件交叉日": sorted(set(cross_days)),
        "口径": "3b:主指标=每层错误率斜率;交叉=裸比特差日的条件性事件,绝不以单日交叉充当结论",
    }


# ---------------------------------------------------------------- e07(三幕)
def e07(cutoff: str) -> dict:
    days = []
    for date, doc in _load("e07_flag", cutoff):
        rows = doc.get("rows", doc) if isinstance(doc, dict) else doc
        if not isinstance(rows, list):
            continue
        rec = {"date": date}
        prov = doc.get("provenance", {}) if isinstance(doc, dict) else {}
        rec["site"] = str(prov.get("selection", prov.get("site", "")))[:60]
        for r in rows:
            if r.get("case") in ("none", "X_anc", "X_D0", "X_D1"):
                rec[r["case"]] = r.get("p_flag_1")
        if "none" in rec and "X_anc" in rec:
            rec["Δ"] = round(rec["X_anc"] - rec["none"], 4)
        days.append(rec)
    acts = {"第一幕_0708_0713": [], "第二幕_0714_0721": [], "第三幕_0722起": []}
    for rec in days:
        d = rec["date"]
        key = ("第一幕_0708_0713" if d <= "2026-07-13" else
               "第二幕_0714_0721" if d <= "2026-07-21" else "第三幕_0722起")
        acts[key].append(rec)
    out = {"逐日": days}
    for act, recs in acts.items():
        out[act] = {
            "天数": len(recs),
            "本底": _pool([r.get("none") for r in recs]),
            "hook": _pool([r.get("X_anc") for r in recs]),
            "Δ": _pool([r.get("Δ") for r in recs]),
        }
    return out


# ---------------------------------------------------------------- e09(树vs线)
def e09(cutoff: str) -> dict:
    per = defaultdict(lambda: defaultdict(list))  # (chip,layout,N)→date→[fid代理]
    for date, doc in _load("e09_ghz", cutoff):
        by = defaultdict(dict)
        for r in doc.get("rows", []):
            by[(r.get("chip", "?"), r["layout"], r["N"])][r["basis"]] = r["value"]
        for key, vals in by.items():
            if "Z" in vals and "X" in vals:
                per[key][date].append((vals["Z"] + vals["X"]) / 2)
    out = defaultdict(dict)
    for (chip, layout, N), day_vals in sorted(per.items()):
        out[chip][f"{layout}_N{N}"] = _daily_then_pool(day_vals)
    for chip in list(out):
        pairs = {}
        for k in out[chip]:
            if k.startswith("tree_"):
                lin = k.replace("tree_", "linear_")
                if lin in out[chip]:
                    t, l = out[chip][k], out[chip][lin]
                    diff = {d: round(t["daily_means"][d] - l["daily_means"][d], 4)
                            for d in t["daily_means"] if d in l["daily_means"]}
                    pairs[k.replace("tree_", "")] = {
                        "tree−linear_逐日": diff,
                        "tree胜日": sum(v > 0 for v in diff.values()),
                        "总日": len(diff)}
        out[chip]["树线对比"] = pairs
    return dict(out)


# ---------------------------------------------------------------- e13(收口素材)
def e13(cutoff: str) -> dict:
    days = []
    for date, doc in _load("e13_mitigation", cutoff):
        rec = {"date": date}
        for r in doc.get("rows", []):
            b = r.get("basis")
            if b in ("ZZ", "XX") and f"raw_{b}" not in rec:
                rec[f"raw_{b}"] = r.get("raw")
                rec[f"twirl_{b}"] = r.get("twirl")
        if "raw_ZZ" in rec and "raw_XX" in rec and None not in (
                rec["raw_ZZ"], rec["raw_XX"], rec.get("twirl_ZZ"), rec.get("twirl_XX")):
            rec["不对称_raw"] = round(abs(rec["raw_XX"] - rec["raw_ZZ"]), 4)
            rec["不对称_twirl"] = round(abs(rec["twirl_XX"] - rec["twirl_ZZ"]), 4)
            rec["twirl压平"] = rec["不对称_twirl"] < rec["不对称_raw"]
        days.append(rec)
    ok = [r for r in days if "twirl压平" in r]
    return {"逐日": days,
            "压平天数": f"{sum(r['twirl压平'] for r in ok)}/{len(ok)}",
            "尖峰日(未压平)": [r["date"] for r in ok if not r["twirl压平"]]}


# ---------------------------------------------------------------- e25 / e26
def e25(cutoff: str) -> dict:
    per: dict = defaultdict(lambda: defaultdict(lambda: {
        "err": defaultdict(list), "acc": defaultdict(list), "bare": defaultdict(list)}))
    for date, doc in _load("e25_iceberg", cutoff):
        for r in doc.get("rows", []):
            chip, k = r.get("chip", "?"), r.get("k")
            if r.get("kind") == "ice" and r.get("err_per_logical") is not None:
                per[chip][k]["err"][date].append(r["err_per_logical"])
                if r.get("acceptance") is not None:
                    per[chip][k]["acc"][date].append(r["acceptance"])
            elif r.get("kind") == "bare" and r.get("err_per_logical") is not None:
                per[chip][k]["bare"][date].append(r["err_per_logical"])
    out = {}
    for chip, by_k in per.items():
        out[chip] = {f"k{k}": {
            "每逻辑比特错误率": _daily_then_pool(v["err"]),
            "接受率": _daily_then_pool(v["acc"]),
            "裸对照": _daily_then_pool(v["bare"]),
            "编码率k/(k+2)": round(k / (k + 2), 3),
        } for k, v in sorted(by_k.items())}
    return out


def e26(cutoff: str) -> dict:
    per: dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
    kinds_seen = set()
    for date, doc in _load("e26_code832", cutoff):
        for r in doc.get("rows", []):
            kind = r.get("kind", "?")
            kinds_seen.add(kind)
            chip = r.get("chip", "?")
            if kind == "ccz" and isinstance(r.get("corr"), (int, float)):
                tag = "corr_withT" if r.get("with_t") else "corr_noT"
                per[chip][kind][tag][date].append(r["corr"] * (r.get("ideal") or 1))
            for field in ("err_any", "acceptance"):
                if isinstance(r.get(field), (int, float)):
                    per[chip][kind][field][date].append(r[field])
    out = {"行类别": sorted(kinds_seen)}
    for chip, by_kind in per.items():
        out[chip] = {kind: {f: _daily_then_pool(dv) for f, dv in fields.items()}
                     for kind, fields in by_kind.items()}
    return out


# ---------------------------------------------------------------- e20 选址A/B
def e20(cutoff: str) -> dict:
    ratio_day = defaultdict(lambda: defaultdict(list))  # d→date→[frozen/daily]
    for date, doc in _load("e20_anchor_ab", cutoff):
        arms = defaultdict(dict)
        for r in doc.get("rows", []):
            if r.get("kind") == "rep" and "err" in r:
                arms[r["d"]][r["arm"]] = r["err"]
        for d, a in arms.items():
            if a.get("daily") and a.get("frozen") is not None and a["daily"] > 1e-6:
                ratio_day[d][date].append(a["frozen"] / a["daily"])
    return {f"d{d}_冻结链err/当日链err": _daily_then_pool(dv)
            for d, dv in sorted(ratio_day.items())}


# ---------------------------------------------------------------- e12 解码阶梯
def e12(cutoff: str) -> dict:
    per = defaultdict(lambda: defaultdict(list))
    for date, doc in _load("e12_decoder", cutoff):
        for r in doc.get("rows", []):
            if r.get("bit") != 1 or r.get("code") != "bit":
                continue
            d = r["d"]
            for field in ("A_majority", "C_mwpm_uniform", "D_mwpm_selfcal",
                          "E_mwpm_empirical", "F_map_q0", "G_map_qfit"):
                if isinstance(r.get(field), (int, float)):
                    tag = field + ("_holdout" if field == "G_map_qfit"
                                   and r.get("G_map_qfit_in_sample") is False else "")
                    per[(d, tag)][date].append(r[field])
    out = defaultdict(dict)
    for (d, field), dv in sorted(per.items()):
        out[f"d{d}"][field] = _daily_then_pool(dv)
    # MWPM 恶化与 MAP 增益(vs 多数表决),取 d 最大档
    dmax = max((int(k[1:]) for k in out), default=None)
    if dmax:
        g = out[f"d{dmax}"]
        if g.get("A_majority", {}).get("n") and g.get("C_mwpm_uniform", {}).get("n"):
            out["判定"] = {
                "d档": dmax,
                "MWPM/表决_错误率比": round(
                    g["C_mwpm_uniform"]["mean"] / g["A_majority"]["mean"], 3)
                if g["A_majority"]["mean"] else None,
                "MAP/表决_错误率比": round(
                    g["F_map_q0"]["mean"] / g["A_majority"]["mean"], 3)
                if g.get("F_map_q0", {}).get("mean") and g["A_majority"]["mean"] else None,
            }
    return dict(out)


# ---------------------------------------------------------------- e06 检错后选
def e06(cutoff: str) -> dict:
    per = defaultdict(lambda: defaultdict(list))
    for date, doc in _load("e06_postselect", cutoff):
        for r in doc.get("rows", []):
            if r.get("bit") != 1 or r.get("code") != "bit":
                continue
            d = r["d"]
            for f in ("raw_majority_err", "ps_err", "acceptance"):
                if isinstance(r.get(f), (int, float)):
                    per[(d, f)][date].append(r[f])
    out = defaultdict(dict)
    for (d, f), dv in sorted(per.items()):
        out[f"d{d}"][f] = _daily_then_pool(dv)
    return dict(out)


# ---------------------------------------------------------------- e17 方差校准
def e17(cutoff: str) -> dict:
    day_ratio = defaultdict(list)
    for date, doc in _load("e17_variance", cutoff):
        rows = doc if isinstance(doc, list) else doc.get("rows", [doc])
        for r in rows:
            if not isinstance(r, dict):
                continue
            for v in r.values():
                er = v.get("excess_ratio") if isinstance(v, dict) else None
                if isinstance(er, (int, float)):
                    day_ratio[date].append(er)
    return {"job间涨落/散粒噪声理论值": _daily_then_pool(day_ratio),
            "注": "≈1 即 job 间涨落与 shot 噪声同量级,误差棒经实测校准"}


# ---------------------------------------------------------------- e18 时刻定位
def e18(cutoff: str) -> dict:
    per = defaultdict(lambda: defaultdict(list))
    hit = defaultdict(lambda: [0, 0])
    for date, doc in _load("e18_t2swap", cutoff):
        for r in doc.get("rows", []):
            case = r.get("case")
            if case and isinstance(r.get("p_expected"), (int, float)):
                per[case][date].append(r["p_expected"])
                dist = r.get("dist") or {}
                if dist:
                    exp_key = "".join(map(str, r.get("expected", [])))
                    top = max(dist, key=dist.get)
                    hit[case][0] += int(top == exp_key or top[::-1] == exp_key)
                    hit[case][1] += 1
    return {case: {"特征概率": _daily_then_pool(dv),
                   "最大概率命中特征": f"{hit[case][0]}/{hit[case][1]}"}
            for case, dv in sorted(per.items())}


# ---------------------------------------------------------------- e08 并行tiling
def e08(cutoff: str) -> dict:
    day_delta = defaultdict(list)
    day_ntiles = defaultdict(list)
    for date, doc in _load("e08_tiling", cutoff):
        prov = doc.get("provenance", {})
        if prov.get("n_tiles"):
            day_ntiles[date].append(prov["n_tiles"])
        for r in doc.get("rows", []):
            if isinstance(r.get("delta"), (int, float)):
                day_delta[date].append(r["delta"])
    return {"packed-isolated_串扰Δ": _daily_then_pool(day_delta),
            "最大同日tile数": max((max(v) for v in day_ntiles.values()), default=None)}


# ---------------------------------------------------------------- e23 / e22 / 漂移
def e23(cutoff: str) -> dict:
    out = defaultdict(lambda: defaultdict(dict))
    for date, doc in _load("e23_tomography", cutoff):
        chip = doc.get("chip") or doc.get("provenance", {}).get("chip", "?")
        for state, blob in (doc.get("states") or {}).items():
            if not isinstance(blob, dict):
                continue
            rec = {k: round(v, 4) for k, v in blob.items()
                   if isinstance(v, (int, float)) and any(
                       t in k.lower() for t in ("fid", "chsh", "mermin", "purity", "neg"))}
            if rec:
                out[chip][state][date] = rec
    return {c: dict(v) for c, v in out.items()}


def e22(cutoff: str) -> dict:
    per = defaultdict(lambda: defaultdict(list))
    for date, doc in _load("e22_variance_panel", cutoff):
        for r in doc.get("rows", []):
            lab = r.get("label") or r.get("kind") or "?"
            for f, v in r.items():
                if isinstance(v, (int, float)) and f not in ("tid", "shots"):
                    per[f"{lab}.{f}"][date].append(v)
    return {k: _daily_then_pool(dv) for k, dv in sorted(per.items()) if len(dv) >= 3}


def drift(cutoff: str, e24_stats: dict) -> dict:
    day_eta = defaultdict(dict)
    pl = R / "e14_probe" / "probe_log.jsonl"
    if pl.exists():
        for line in pl.open():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            d = rec.get("date", "")
            if FIRST_DAY <= d <= cutoff and rec.get("selfcal", {}).get("eta_median"):
                day_eta[rec.get("chip", "?")][d] = rec["selfcal"]["eta_median"]
    out = {"eta_median_逐日": {c: dict(sorted(v.items())) for c, v in day_eta.items()}}
    bh = day_eta.get("Baihua", {})
    dm = (e24_stats.get("Baihua") or {}).get("跨日_逻辑F", {}).get("daily_means", {})
    common = sorted(set(bh) & set(dm))
    if len(common) >= 3:
        out["漂移-性能相关(附录B 2)"] = {
            "eta_median vs e24逻辑F 日级r": _pearson(
                [bh[d] for d in common], [dm[d] for d in common]),
            "n_days": len(common),
            "注": "机器状态可预测实验质量;轮级证据=e24 接受率哨兵",
        }
    return out


# ---------------------------------------------------------------- md 渲染
def _fmt(p: dict) -> str:
    if not p or not p.get("n"):
        return "—"
    if p["n"] == 1:
        return f"{p['mean']:.4f}(n=1)"
    return f"{p['mean']:.4f} [{p['lo']:.4f},{p['hi']:.4f}] (n={p['n']}, 中位 {p['median']:.4f})"


def render_md(S: dict) -> str:
    L = [f"# 终局统计数字表(统计截至 {S['meta']['cutoff']},生成于 {S['meta']['generated']})",
         "", "> 本表由 tools/freeze_stats.py 自动生成;正文引用数字一律以本表为准。",
         "> 口径:跨日均值 ± bootstrap 95% CI,n=天数;后选必附接受率。", ""]
    w = S["workload"]
    L += ["## 真机工作量(评审维度:真机实践)", "",
          f"- 有效天数 **{w['有效天数']}**({w['起止'][0]} → {w['起止'][1]}),三机:{'/'.join(w['机器'])}",
          f"- 提交任务总数 **{w['总任务数']}**(日均 {w['日均任务']});summary 文件 {w['summary文件总数']} 个", ""]
    L += ["## e24 逻辑 Bell vs 物理 Bell(headline #1)", ""]
    for chip, s in S["e24"].items():
        L += [f"### {chip}({s['轮数']} 轮 / {s['天数']} 天,日胜负 {s['日胜负']})", "",
              f"- 跨日逻辑 F:{_fmt(s['跨日_逻辑F'])}",
              f"- 跨日物理 F:{_fmt(s['跨日_物理F'])}",
              f"- 跨日接受率:{_fmt(s['跨日_接受率'])}",
              f"- 哨兵相关:全轮 r={s['哨兵_全轮r']},日级 r={s['哨兵_日级r']}"]
        st = s.get("双态分层")
        if st:
            L += [f"- **双态分层(阈值 {st['阈值(2means)']},坏态占比 {st['坏态占比']:.0%})**:",
                  f"  - 好态({st['好态轮数']} 轮):逻辑 {_fmt(st['好态_逻辑F'])} vs 物理 {_fmt(st['好态_物理F'])}",
                  f"  - 坏态({st['坏态轮数']} 轮):逻辑 {_fmt(st['坏态_逻辑F'])} vs 物理 {_fmt(st['坏态_物理F'])}"]
        L.append("")
    L += ["## e05 距离扫描(headline #2)", ""]
    for chip, curves in S["e05"].items():
        for code in ("bit", "phase"):
            if code in curves:
                row = " | ".join(f"{d}:{_fmt(p).split(' [')[0]}"
                                 for d, p in curves[code].items())
                L.append(f"- {chip} {code}:{row}")
        fl = curves.get("地板判定")
        if fl:
            L.append(f"  - 地板:d9/d7={fl['d9/d7']}({fl['说明']})")
    L.append("")
    L += ["## e11 编码存储对决(headline #3)", ""]
    for chip, s in S["e11"].items():
        if chip == "换址日期":
            L.append(f"- 换址日期记录(3g 口径预答):{s}")
            continue
        for code in ("bit", "phase"):
            hd = s.get(f"headline_d7对决_{code}")
            if hd:
                L += [f"- {chip} {code} d7 {hd['K档']}:enc {_fmt(hd['enc'])} vs "
                      f"最好物理臂 {_fmt(hd['phys_best'])}(CI 分离:{hd['CI分离']})"]
    L.append("")
    L += ["## e10 [[4,2,2]] 对决 + 融合制备", "",
          f"- 最好星(后选 perq):{_fmt(S['e10']['最好星_后选perq'])},接受率 {_fmt(S['e10']['最好星_接受率'])}",
          f"- 最强物理基线:{_fmt(S['e10']['最强物理基线_perq'])}",
          f"- 融合制备:{_fmt(S['e10']['融合制备_perq'])},接受率 {_fmt(S['e10']['融合_接受率'])}",
          f"- 判定:{S['e10'].get('对决判定')};融合:{S['e10'].get('融合判定')}", ""]
    e19s = S["e19"]
    L += ["## e19 门级斜率(3b 口径)", "",
          f"- enc 斜率:{_fmt(e19s['enc斜率_跨日'])}",
          f"- bare 斜率:{_fmt(e19s['bare斜率_跨日'])};enc 更低天数 {e19s['enc斜率更低的天数']}",
          f"- 斜率比(引用中位数):{_fmt(e19s['斜率比_跨日'])}",
          f"- 条件交叉日:{e19s['条件交叉日']}", ""]
    L += ["## e07 flag 三幕", ""]
    for act in ("第一幕_0708_0713", "第二幕_0714_0721", "第三幕_0722起"):
        a = S["e07"][act]
        L.append(f"- {act}({a['天数']} 天):本底 {_fmt(a['本底'])} / hook {_fmt(a['hook'])} / Δ {_fmt(a['Δ'])}")
    L.append("")
    L += ["## e13 twirl 收口(3g)", "",
          f"- 不对称压平天数:{S['e13']['压平天数']};尖峰日:{S['e13']['尖峰日(未压平)']}", ""]
    L += ["## e09 GHZ 树vs线(3g)", ""]
    for chip, s in S["e09"].items():
        cmp_ = s.get("树线对比", {})
        for N, v in cmp_.items():
            L.append(f"- {chip} {N}:tree 胜 {v['tree胜日']}/{v['总日']} 日")
    L.append("")
    L += ["## e20 自适应选址增益(同批 A/B)", ""]
    for k, p in S["e20"].items():
        L.append(f"- {k}:{_fmt(p)}")
    L.append("")
    j = S["e12"].get("判定", {})
    L += ["## e12 解码阶梯", "",
          f"- 判定(d{j.get('d档')}):MWPM/表决 错误率比 {j.get('MWPM/表决_错误率比')};"
          f"MAP/表决 {j.get('MAP/表决_错误率比')}", ""]
    L += ["## e06 检错后选", ""]
    for dk in ("d3", "d7"):
        g = S["e06"].get(dk)
        if g:
            L.append(f"- {dk}:raw {_fmt(g.get('raw_majority_err', {}))} → ps "
                     f"{_fmt(g.get('ps_err', {}))},接受率 {_fmt(g.get('acceptance', {}))}")
    L.append("")
    L += ["## e17 方差校准", "", f"- {_fmt(S['e17']['job间涨落/散粒噪声理论值'])}({S['e17']['注']})", ""]
    L += ["## e18 两轮时刻定位", ""]
    for case, g in S["e18"].items():
        L.append(f"- {case}:特征概率 {_fmt(g['特征概率'])},命中 {g['最大概率命中特征']}")
    L.append("")
    L += ["## 漂移-性能(附录 B 2)", "", f"- {S['drift'].get('漂移-性能相关(附录B 2)', '—')}", ""]
    if WARN:
        L += ["## ⚠️ 引擎警告", ""] + [f"- {wn}" for wn in WARN] + [""]
    L += ["", "*(e22/e23/e25/e26 全量数字在 stats.json 对应节;本表只列报告要引用的主线。)*"]
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cutoff", default=None, help="统计截止日 YYYY-MM-DD,默认今天")
    args = ap.parse_args()
    import datetime
    cutoff = args.cutoff or datetime.date.today().isoformat()
    gen = datetime.datetime.now().strftime("%F %H:%M")

    S = {"meta": {"cutoff": cutoff, "generated": gen, "first_day": FIRST_DAY}}
    S["workload"] = workload(cutoff)
    S["e24"] = e24(cutoff)
    S["e05"] = e05(cutoff)
    S["e11"] = e11(cutoff)
    S["e10"] = e10(cutoff)
    S["e19"] = e19(cutoff)
    S["e07"] = e07(cutoff)
    S["e09"] = e09(cutoff)
    S["e13"] = e13(cutoff)
    S["e25"] = e25(cutoff)
    S["e26"] = e26(cutoff)
    S["e08"] = e08(cutoff)
    S["e20"] = e20(cutoff)
    S["e12"] = e12(cutoff)
    S["e06"] = e06(cutoff)
    S["e17"] = e17(cutoff)
    S["e18"] = e18(cutoff)
    S["e23"] = e23(cutoff)
    S["e22"] = e22(cutoff)
    S["drift"] = drift(cutoff, S["e24"])
    S["warnings"] = WARN

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "stats.json").write_text(json.dumps(S, ensure_ascii=False, indent=1))
    (OUT / "stats.md").write_text(render_md(S))
    print(f"已写 {OUT}/stats.json + stats.md(截至 {cutoff})")
    for wn in WARN:
        print(f"  ⚠️ {wn}")


if __name__ == "__main__":
    main()
