"""e15 · 跨机 capability map(Baihua 榨干后启动;后端上线窗口机会主义运行)。

对指定后端跑迷你套件:全芯片自校准(2)+ 自动选链 repcode d=3 L0/L1(2)
+ 最优边 Bell ZZ/XX(2)= 6 任务,产出与 Baihua 同口径的对比行。

运行: ./qrun experiments/e15_crossbackend.py <chip>   # 如 Dongling / Shenglian
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, chip as chipmod, client, provenance, repcode  # noqa: E402
from qfab.qasm import QASM, bell_pair  # noqa: E402

TODAY = time.strftime("%Y-%m-%d")
PROJ = Path(__file__).resolve().parent.parent
OUT = PROJ / "results" / "e15_crossbackend"


def auto_chain(ch, n: int = 5):
    """按 vendor 标定贪心找一条 n 比特链(自测数据要等 selfcal 之后)。"""
    av = set(ch.available_qubits())
    def score(q):
        return ch.qubits[q].get("readout_fidelity") or 0
    best = None
    for start in av:
        path = [start]
        while len(path) < n:
            cands = [x for x in ch.neighbors(path[-1]) if x in av
                     and x not in path and (ch.edge_fid(path[-1], x) or 0) > 0.5]
            if not cands:
                break
            path.append(max(cands, key=score))
        if len(path) == n:
            s = min(score(q) for q in path)
            if best is None or s > best[0]:
                best = (s, path)
    return best[1] if best else None


def main():
    if len(sys.argv) < 2:
        print("用法: ./qrun experiments/e15_crossbackend.py <chip>")
        return
    chip = sys.argv[1]
    st = client.status()
    print(f"后端状态: {st}")
    if not isinstance(st.get(chip), int):
        print(f"⚠️ {chip} 不可用({st.get(chip)})")
        return
    ch = chipmod.load(chip, refresh=True)
    print(f"芯片: {ch.summary()}")
    av = ch.available_qubits()
    chain = auto_chain(ch)
    edges = sorted(((f or 0, a, b) for a, b, f in ch.edges
                    if a in set(av) and b in set(av)), reverse=True)
    best_edge = (edges[0][1], edges[0][2]) if edges else None
    print(f"自动链: {chain},最优边: {best_edge}")

    exp = f"e15_crossbackend/{chip}_{TODAY}"
    # Shenglian 单任务宽度墙(>32~40 比特挂死)→ 自校准切块;
    # 且只测厂商读出已标定(g/e>0)的比特,读出未标定的测了也是噪声还可能挂任务
    max_w = {"Shenglian": 28}.get(chip, 10**9)
    sc = av
    if len(sc) > max_w:
        def _rd(q):
            raw = ch.qubits[q]["raw"]
            return ((raw.get("readout g_fidelity") or 0)
                    + (raw.get("readout e_fidelity") or 0))
        sc = [q for q in sc if _rd(q) > 0]
    groups = [sc[i:i + max_w] for i in range(0, len(sc), max_w)]
    specs = []
    for gi, g in enumerate(groups):
        n = len(g)
        qc0 = QASM(n); qc0.measure_all()
        qc1 = QASM(n)
        for i in range(n):
            qc1.x(i)
        qc1.measure_all()
        specs += [
            dict(circuit=qc0.qasm(), chip=chip, shots=4096, name=f"xsc0_g{gi}",
                 exp=exp, label=f"selfcal_prep0_g{gi}",
                 options={"target_qubits": g},
                 meta={"kind": "selfcal", "prep": "prep0", "qubits": g}),
            dict(circuit=qc1.qasm(), chip=chip, shots=4096, name=f"xsc1_g{gi}",
                 exp=exp, label=f"selfcal_prep1_g{gi}",
                 options={"target_qubits": g},
                 meta={"kind": "selfcal", "prep": "prep1", "qubits": g}),
        ]
    if chain:
        for bit in (0, 1):
            specs.append(dict(circuit=repcode.build_t1(3, bit), chip=chip,
                              shots=4096, name=f"xrep{bit}", exp=exp,
                              label=f"rep_L{bit}",
                              options={"target_qubits": chain},
                              meta={"kind": "repcode", "d": 3, "bit": bit,
                                    "chain": chain}))
    if best_edge:
        for basis in ("ZZ", "XX"):
            specs.append(dict(circuit=bell_pair(basis).qasm(), chip=chip,
                              shots=2048, name=f"xbell{basis}", exp=exp,
                              label=f"bell_{basis}",
                              options={"target_qubits": list(best_edge)},
                              meta={"kind": "bell", "edge": list(best_edge),
                                    "basis": basis}))
    recs = client.run_batch(specs, timeout=1200)

    by = {r["label"]: r for r in recs if client.counts(r)}
    day = {"date": TODAY, "chip": chip, "chain": chain, "edge": best_edge}
    etas = []
    for gi, g in enumerate(groups):
        c0, c1 = (client.counts(by.get(f"selfcal_prep{p}_g{gi}", {}))
                  for p in (0, 1))
        if c0 and c1:
            asg = analysis.assignment_from_p1(analysis.bit_p1(c0, len(g)),
                                              analysis.bit_p1(c1, len(g)))
            etas += [a["eta"] for a in asg]
    if etas:
        etas.sort()
        day["eta_median"] = round(etas[len(etas) // 2], 4)
        day["n_selfcal_qubits"] = len(etas)
    for bit in (0, 1):
        r = by.get(f"rep_L{bit}")
        if r:
            c = client.counts(r)
            nn = analysis.shots_of(c)
            me = repcode.logical_error_majority(c, 3, bit)
            lo, hi = analysis.wilson_ci(round(me * nn), nn)
            day[f"rep_d3_L{bit}"] = round(me, 5)
            day[f"rep_d3_L{bit}_ci"] = [round(lo, 5), round(hi, 5)]
    for basis in ("ZZ", "XX"):
        r = by.get(f"bell_{basis}")
        if r:
            c = client.counts(r)
            nn = analysis.shots_of(c)
            even = sum(v for k, v in c.items() if k.zfill(2).count("1") % 2 == 0)
            day[f"bell_{basis}"] = round(2 * even / nn - 1, 4)
            lo, hi = analysis.wilson_ci(even, nn)
            day[f"bell_{basis}_ci"] = [round(2 * lo - 1, 4), round(2 * hi - 1, 4)]
    # 溯源戳:注意选址口径=vendor 读出保真(auto_chain/最优 CZ 边),与主线 Baihua 自测 η 口径不同。
    day["provenance"] = provenance.stamp(
        anchor=None, qubits=sorted(set(chain or []) | set(best_edge or ())),
        tids=provenance.tids_of(recs), chip=chip,
        extra={"chip_calibration_time": ch.calibration_time,
               "site_selection": "vendor_readout_fidelity (auto_chain / best-CZ edge)",
               "caliber_note": "选址口径=vendor 标定读出保真,与主线 Baihua 自测 η 口径不一致(跨机榨干前无自测数据)"})
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "crossbackend_log.jsonl").open("a", encoding="utf-8") as f:
        f.write(json.dumps(day, ensure_ascii=False) + "\n")
    print(json.dumps(day, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
