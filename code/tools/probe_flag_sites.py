"""e07 flag 选址探针(一次性,10 任务)。

现状:e07 硬编码 syndrome 锚(data[20,22]/anc[21]),flag 取 anc 第三邻居
(edge>0.5 即收)——该环节从未被选址,flag 本底 P(flag=1|none) 高达 0.4-0.65,
hook 信号被完全淹没。本探针全芯片枚举候选星(anc 有 ≥3 可用邻居),按
flag 环节质量 edge_fid(anc,flag)² × readout(flag) 排序,取 top4 + 现任站
(对照),每站只测判别对 none / X_anc(2 任务 × 5 站 = 10 任务,探针纪律内)。

判据:Δ = P(flag=1 | X_anc) − P(flag=1 | none)。
可用站 = 本底 < 0.15 且 Δ 的 Wilson 区间显著 > 0;若全芯片无可用站,
即为量化负结果:flag 原语入场门槛高于此机当前边质量。

运行: ./qrun tools/probe_flag_sites.py [chip]     # 默认 Baihua
      ./qrun tools/probe_flag_sites.py selftest    # 零额度 statevector 自检
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))

from experiments.e07_flag import circ  # noqa: E402  复用 e07 真实线路,不另写物理
from qfab import analysis, chip as chipmod, client, provenance  # noqa: E402

CHIP = "Baihua"
RUN = time.strftime("%Y-%m-%d_%H%M")
OUT = PROJ / "results" / "e07_flag"
SHOTS = 4096
CASES = ("none", "X_anc")            # 判别对;X_D0/X_D1 留给定址后的 e07 正式轮
TOPN = 4


def selftest() -> int:
    from tools.verify_circuits_offline import simulate
    exp = {"none": (0, 0), "X_anc": (1, 1), "X_D0": (0, 1), "X_D1": (0, 1)}
    ok = True
    for case, (want_flag, want_anc) in exp.items():
        dist = simulate(circ(case))
        p_flag = sum(p for c, p in dist.items() if c[0] == "1")   # c3=flag(最左)
        p_anc = sum(p for c, p in dist.items() if c[1] == "1")    # c2=anc
        good = abs(p_flag - want_flag) < 1e-9 and abs(p_anc - want_anc) < 1e-9
        ok &= good
        print(f"  {case:6s}: P(flag=1)={p_flag:.0f} P(anc=1)={p_anc:.0f} "
              f"期望=({want_flag},{want_anc}) {'✓' if good else '✗'}")
    print("selftest:", "全绿" if ok else "失败")
    return 0 if ok else 1


def enumerate_sites(ch) -> list[dict]:
    av = set(ch.available_qubits())

    def ro(q):
        return (ch.qubits.get(q) or {}).get("readout_fidelity") or 0

    def flag_score(anc, x):
        return (ch.edge_fid(anc, x) or 0) ** 2 * ro(x)

    sites = []
    for anc in sorted(av):
        nb = [x for x in ch.neighbors(anc) if x in av and (ch.edge_fid(anc, x) or 0) > 0]
        if len(nb) < 3:
            continue
        fl = max(nb, key=lambda x: flag_score(anc, x))
        rest = sorted((x for x in nb if x != fl),
                      key=lambda x: ch.edge_fid(anc, x) or 0, reverse=True)
        sites.append({"anc": anc, "flag": fl, "data": rest[:2],
                      "flag_edge": round(ch.edge_fid(anc, fl), 4),
                      "flag_readout": round(ro(fl), 4),
                      "score": round(flag_score(anc, fl), 4)})
    sites.sort(key=lambda s: s["score"], reverse=True)
    return sites


def incumbent_site(ch) -> dict | None:
    """复刻 e07 现行规则(data[20,22]/anc21,第三邻居 edge>0.5 取首个)做对照。"""
    av = set(ch.available_qubits())
    flags = [x for x in ch.neighbors(21)
             if x in av and x not in (20, 22) and (ch.edge_fid(21, x) or 0) > 0.5]
    if not flags or not {20, 21, 22} <= av:
        return None
    fl = flags[0]
    ro = (ch.qubits.get(fl) or {}).get("readout_fidelity") or 0
    return {"anc": 21, "flag": fl, "data": [20, 22],
            "flag_edge": round(ch.edge_fid(21, fl) or 0, 4),
            "flag_readout": round(ro, 4),
            "score": round((ch.edge_fid(21, fl) or 0) ** 2 * ro, 4),
            "incumbent": True}


def main():
    ch = chipmod.load(CHIP, refresh=True)
    sites = enumerate_sites(ch)[:TOPN]
    inc = incumbent_site(ch)
    if inc and not any(s["anc"] == inc["anc"] and s["flag"] == inc["flag"] for s in sites):
        sites.append(inc)
    if not sites:
        print("⚠️ 无候选站(可用邻居不足),中止")
        return
    print(f"候选站 ×{len(sites)}(按 flag 环节质量降序,末位带 incumbent=现任对照):")
    for s in sites:
        print(f"  anc=Q{s['anc']:<3} flag=Q{s['flag']:<3} data={s['data']} "
              f"边={s['flag_edge']} 读出={s['flag_readout']} 分={s['score']}"
              + (" [现任]" if s.get("incumbent") else ""))

    exp = f"e07_flag/probe_{RUN}"
    specs = []
    for s in sites:
        tq = [*s["data"], s["anc"], s["flag"]]
        tag = f"a{s['anc']}f{s['flag']}"
        for cs in CASES:
            specs.append(dict(circuit=circ(cs), chip=CHIP, shots=SHOTS,
                              name=f"probe_{tag}_{cs}", exp=exp, label=f"{tag}_{cs}",
                              options={"target_qubits": tq},
                              meta={"site": tag, "case": cs, **{k: s[k] for k in
                                    ("anc", "flag", "data", "flag_edge",
                                     "flag_readout", "score")}}))
    recs = client.run_batch(specs, timeout=1200)

    by_site: dict[str, dict] = {}
    for r in recs:
        c = client.counts(r)
        if not c:
            continue
        m = r["meta"]
        p1 = analysis.bit_p1(c, 4)
        tot = analysis.shots_of(c)
        ci = analysis.wilson_ci(round(p1[3] * tot), tot)
        by_site.setdefault(m["site"], {k: m[k] for k in
                           ("anc", "flag", "data", "flag_edge", "flag_readout")})
        by_site[m["site"]][m["case"]] = {
            "p_flag_1": round(p1[3], 4),
            "p_flag_1_ci": [round(ci[0], 4), round(ci[1], 4)],
            "p_anc_1": round(p1[2], 4), "shots": tot}

    rows = []
    print("\n===== 判别力(Δ = P(flag=1|X_anc) − P(flag=1|none))=====")
    for tag, d in by_site.items():
        if "none" not in d or "X_anc" not in d:
            print(f"  {tag}: 结果不全(队列超时?先 recollect 再看)")
            rows.append({"site": tag, **d, "verdict": "incomplete"})
            continue
        base, sig = d["none"], d["X_anc"]
        delta = round(sig["p_flag_1"] - base["p_flag_1"], 4)
        usable = base["p_flag_1"] < 0.15 and sig["p_flag_1_ci"][0] > base["p_flag_1_ci"][1]
        verdict = "usable" if usable else "noisy"
        rows.append({"site": tag, **d, "delta": delta, "verdict": verdict})
        print(f"  {tag}(边{d['flag_edge']}/读出{d['flag_readout']}): "
              f"本底={base['p_flag_1']:.4f} hook={sig['p_flag_1']:.4f} "
              f"Δ={delta:+.4f} anc核对={sig['p_anc_1']:.4f} → {verdict}")

    good = [r for r in rows if r["verdict"] == "usable"]
    if good:
        best = max(good, key=lambda r: r["delta"])
        print(f"\n推荐站: {best['site']}(Δ={best['delta']:+.4f});"
              f"改 e07 选址后跑 selftest + verify_circuits_offline 再上真机")
    else:
        print("\n无可用站 → 负结果素材:flag 原语入场门槛高于此机当前边质量")

    prov = provenance.stamp(
        qubits=sorted({q for s in sites for q in (*s["data"], s["anc"], s["flag"])}),
        tids=provenance.tids_of(recs), chip=CHIP,
        extra={"selection": "flag_score = edge(anc,flag)^2 × readout(flag), top4+incumbent",
               "cases": list(CASES), "shots_per_task": SHOTS})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"probe_sites_{RUN}.json").write_text(
        json.dumps({"provenance": prov, "rows": rows}, ensure_ascii=False, indent=2))
    print(f"已写 {OUT / f'probe_sites_{RUN}.json'}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "selftest":
        sys.exit(selftest())
    CHIP = sys.argv[1] if len(sys.argv) > 1 else CHIP
    main()
