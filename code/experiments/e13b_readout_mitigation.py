"""e13b · 读出误差缓解演示(标准方法,纯离线,零额度)。

用当日自测 assignment 矩阵(e02)对已有真机数据做张量化反演缓解,三个演示对象:
  A. repetition 距离扫描(e05):raw vs 缓解后的多数表决逻辑错误率
  B. GHZ 布居(e09):raw vs 缓解后的 P(全0)+P(全1)
  C. [[4,2,2]](e10 最好星):raw/后选/缓解/缓解+后选 四档对照
均报告负值截断占比(反演病态度的诚实指标)。

运行: ./qrun experiments/e13b_readout_mitigation.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, mitigation, provenance, repcode  # noqa: E402

PROJ = Path(__file__).resolve().parent.parent
OUT = PROJ / "results" / "e13_mitigation"
RUN = time.strftime("%Y-%m-%d_%H%M")


def load_cal() -> dict[int, dict]:
    files = sorted((PROJ / "results" / "e02_selfcal").glob("assignment_*.json"))
    return {r["qubit"]: r for r in json.loads(files[-1].read_text())}


def latest_dir(exp: str) -> Path:
    return sorted(d for d in (PROJ / "results" / exp).iterdir() if d.is_dir())[-1]


def main():
    cal = load_cal()
    cal_src = sorted((PROJ / "results" / "e02_selfcal").glob("assignment_*.json"))[-1].name
    rows = []
    used_recs = []
    d05, d09, d10 = (latest_dir("e05_repcode"), latest_dir("e09_ghz"),
                     latest_dir("e10_code422"))

    print("===== A. repetition:raw vs 读出缓解(多数表决逻辑错误率)=====")
    for f in sorted(d05.glob("bit_d*_L*__*.json")):
        rec = json.loads(f.read_text())
        c = (rec.get("result") or {}).get("count")
        if not c:
            continue
        used_recs.append(rec)
        m = rec["meta"]
        tot = analysis.shots_of(c)
        fs = mitigation.fs_for_qubits(rec["task"]["options"]["target_qubits"], cal)
        cm, clip = mitigation.correct_counts(c, fs)
        raw = repcode.logical_error_majority(c, m["d"], m["bit"])
        mit = repcode.logical_error_majority(cm, m["d"], m["bit"])
        gain = raw / mit if mit > 0 else float("inf")
        raw_lo, raw_hi = analysis.wilson_ci(round(raw * tot), tot)   # 真实计数,精确
        mit_lo, mit_hi = analysis.wilson_ci(round(mit * tot), tot)   # 伪计数,名义
        rows.append({"part": "A", **{k: m[k] for k in ("d", "bit")},
                     "raw": round(raw, 5), "raw_ci": [round(raw_lo, 5), round(raw_hi, 5)],
                     "mitigated": round(mit, 5),
                     "mit_ci_nominal": [round(mit_lo, 5), round(mit_hi, 5)],
                     "gain": round(gain, 2), "clipped": round(clip, 4)})
        print(f"  d={m['d']} L{m['bit']}: {raw:.5f} → {mit:.5f} "
              f"({gain:.1f}×, 截断 {clip:.3f})")

    print("\n===== B. GHZ 布居:raw vs 读出缓解 =====")
    for f in sorted(d09.glob("tree_N*_Z__*.json")):
        rec = json.loads(f.read_text())
        c = (rec.get("result") or {}).get("count")
        if not c:
            continue
        used_recs.append(rec)
        n = rec["meta"]["N"]
        tot = analysis.shots_of(c)
        fs = mitigation.fs_for_qubits(rec["meta"]["qubits"], cal)
        cm, clip = mitigation.correct_counts(c, fs)

        def pops(cc):
            t = analysis.shots_of(cc) or 1
            return sum(v for k, v in cc.items()
                       if set(k.zfill(n)) <= {"0"} or set(k.zfill(n)) <= {"1"}) / t
        p_raw, p_mit = pops(c), pops(cm)
        raw_lo, raw_hi = analysis.wilson_ci(round(p_raw * tot), tot)
        mit_lo, mit_hi = analysis.wilson_ci(round(p_mit * tot), tot)  # 伪计数,名义
        rows.append({"part": "B", "N": n, "raw": round(p_raw, 4),
                     "raw_ci": [round(raw_lo, 4), round(raw_hi, 4)],
                     "mitigated": round(p_mit, 4),
                     "mit_ci_nominal": [round(mit_lo, 4), round(mit_hi, 4)],
                     "clipped": round(clip, 4)})
        print(f"  树形 N={n}: pop {p_raw:.4f} → {p_mit:.4f} (截断 {clip:.3f})")

    print("\n===== C. [[4,2,2]] 四档对照(最好星,Z 基 4 态平均)=====")
    recs = [json.loads(f.read_text()) for f in d10.glob("c*_Z_*__*.json")]
    recs = [r for r in recs if (r.get("result") or {}).get("count")]
    used_recs += recs
    by_star: dict[int, list] = {}
    for r in recs:
        by_star.setdefault(r["meta"]["star"]["center"], []).append(r)

    def eval4(cc, m):
        tot = acc = e_ps = e_raw = 0
        for key, v in cc.items():
            b = analysis.bits_of(key, 4)
            l1, l2 = b[0] ^ b[2], b[0] ^ b[1]
            tot += v
            e_raw += v * ((l1 != m["m"]) + (l2 != m["n"])) / 2
            if (b[0] ^ b[1] ^ b[2] ^ b[3]) == 0:
                acc += v
                e_ps += v * ((l1 != m["m"]) + (l2 != m["n"])) / 2
        return e_raw / (tot or 1), e_ps / (acc or 1), acc / (tot or 1)

    def counts4(cc, m):
        """池化 CI 用计数:逻辑位错误数(=2×e_raw)、接受数、总 shots。"""
        tot = acc = e_raw2 = e_ps2 = 0.0
        for key, v in cc.items():
            b = analysis.bits_of(key, 4)
            l1, l2 = b[0] ^ b[2], b[0] ^ b[1]
            tot += v
            e_raw2 += v * ((l1 != m["m"]) + (l2 != m["n"]))
            if (b[0] ^ b[1] ^ b[2] ^ b[3]) == 0:
                acc += v
                e_ps2 += v * ((l1 != m["m"]) + (l2 != m["n"]))
        return tot, acc, e_raw2, e_ps2

    best = None
    for center, rs in by_star.items():
        vals = {"raw": 0.0, "ps": 0.0, "mit": 0.0, "mit_ps": 0.0, "acc": 0.0, "clip": 0.0}
        pool = {"tot": 0.0, "acc": 0.0, "er2": 0.0, "ep2": 0.0,
                "mtot": 0.0, "macc": 0.0, "mer2": 0.0, "mep2": 0.0}
        for r in rs:
            m = r["meta"]
            c = r["result"]["count"]
            fs = mitigation.fs_for_qubits([m["star"]["center"], *m["star"]["leaves"]], cal)
            cm, clip = mitigation.correct_counts(c, fs)
            r_raw, r_ps, r_acc = eval4(c, m)
            m_raw, m_ps, _ = eval4(cm, m)
            vals["raw"] += r_raw / len(rs); vals["ps"] += r_ps / len(rs)
            vals["mit"] += m_raw / len(rs); vals["mit_ps"] += m_ps / len(rs)
            vals["acc"] += r_acc / len(rs); vals["clip"] += clip / len(rs)
            t, a, er2, ep2 = counts4(c, m)
            pool["tot"] += t; pool["acc"] += a; pool["er2"] += er2; pool["ep2"] += ep2
            mt, ma, mer2, mep2 = counts4(cm, m)
            pool["mtot"] += mt; pool["macc"] += ma; pool["mer2"] += mer2; pool["mep2"] += mep2
        vals["_pool"] = pool
        if best is None or vals["ps"] < best[1]["ps"]:
            best = (center, vals)
    c0, v = best
    pool = v.pop("_pool")

    def wl(k, n, nd=5):
        lo, hi = analysis.wilson_ci(round(k), round(n))
        return [round(lo, nd), round(hi, nd)]
    rows.append({"part": "C", "star": c0,
                 **{k: round(x, 5) for k, x in v.items()},
                 "acc_ci": wl(pool["acc"], pool["tot"]),
                 "raw_ci": wl(pool["er2"], 2 * pool["tot"]),
                 "ps_ci": wl(pool["ep2"], 2 * pool["acc"]),
                 "mit_ci_nominal": wl(pool["mer2"], 2 * pool["mtot"]),
                 "mit_ps_ci_nominal": wl(pool["mep2"], 2 * pool["macc"]),
                 "ci_method": "pooled counts (逻辑位 n=2×shots);nominal=缓解伪计数"})
    print(f"  star c{c0}: raw={v['raw']:.4f} | 仅后选={v['ps']:.4f} (acc {v['acc']:.2f}) | "
          f"仅缓解={v['mit']:.4f} | 缓解+后选={v['mit_ps']:.4f} (截断 {v['clip']:.3f})")

    print("\n===== D. 公平对照:缓解后的物理基线(与 C 的缓解+后选同权对比)=====")
    phys_recs = [json.loads(f.read_text()) for f in d10.glob("phys_q*__*.json")]
    per_q: dict[int, list] = {}
    for r in phys_recs:
        c = (r.get("result") or {}).get("count")
        if not c:
            continue
        per_q.setdefault(r["meta"]["qubit"], []).append(r)
    used_recs += [r for rs in per_q.values() for r in rs]
    best_mit_phys = None
    for q, rs in per_q.items():
        errs_raw, errs_mit = [], []
        k_raw = n_raw = 0
        k_mit = n_mit = 0.0
        for r in rs:
            c = r["result"]["count"]
            tot = analysis.shots_of(c)
            fs = mitigation.fs_for_qubits([q], cal)
            cm, _ = mitigation.correct_counts(c, fs)
            p1 = analysis.bit_p1(c, 1)[0]
            e_raw = p1 if r["meta"]["bit"] == 0 else 1 - p1
            errs_raw.append(e_raw)
            k_raw += round(e_raw * tot); n_raw += tot
            p1m = analysis.bit_p1(cm, 1)[0]
            e_mit = p1m if r["meta"]["bit"] == 0 else 1 - p1m
            errs_mit.append(e_mit)
            mtot = analysis.shots_of(cm)
            k_mit += e_mit * mtot; n_mit += mtot
        raw_m = sum(errs_raw) / len(errs_raw)
        mit_m = sum(errs_mit) / len(errs_mit)
        raw_lo, raw_hi = analysis.wilson_ci(k_raw, n_raw)
        mit_lo, mit_hi = analysis.wilson_ci(round(k_mit), round(n_mit))  # 伪计数,名义
        print(f"  q{q}: raw={raw_m:.4f} → 缓解后={mit_m:.4f}")
        rows.append({"part": "D", "qubit": q, "raw": round(raw_m, 5),
                     "raw_ci": [round(raw_lo, 5), round(raw_hi, 5)],
                     "mitigated": round(mit_m, 5),
                     "mit_ci_nominal": [round(mit_lo, 5), round(mit_hi, 5)]})
        if best_mit_phys is None or mit_m < best_mit_phys[1]:
            best_mit_phys = (q, mit_m)
    if best_mit_phys:
        print(f"  → 最严苛缓解后物理基线 = q{best_mit_phys[0]}: {best_mit_phys[1]:.4f}"
              f"(与 C 档'缓解+后选'同权对比)")

    chip = used_recs[0].get("chip", "Baihua") if used_recs else "Baihua"
    prov = provenance.stamp(
        chip=chip, tids=provenance.tids_of(used_recs),
        extra={"kind": "offline_reanalysis", "cal_source": cal_src,
               "source_dirs": {"e05_repcode": d05.name, "e09_ghz": d09.name,
                               "e10_code422": d10.name},
               "n_records": len(used_recs),
               "ci": "raw 类=Wilson 95% 精确(真实计数);*_ci_nominal=缓解伪计数 shot-noise "
                     "名义下界(反演放大真实方差);C 档为 4 态池化计数"})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"readout_mit_{RUN}.json").write_text(
        json.dumps({"provenance": prov, "rows": rows}, ensure_ascii=False, indent=2))
    print(f"\n已写 {OUT / f'readout_mit_{RUN}.json'}")


if __name__ == "__main__":
    main()
