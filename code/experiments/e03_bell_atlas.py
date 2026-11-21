"""e03 · 全芯片 Bell 边图谱(边着色并行:每个颜色类的所有不相交边塞进一条线路)。

每条标定边测 Bell 态 ZZ/XX parity → 边质量自测图谱 + vendor CZ 对照 +
ZZ/XX 不对称(CZ 相干相位误差指纹,e13 twirling 的对照基线)。

运行: ./qrun experiments/e03_bell_atlas.py [chip]
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, chip as chipmod, client, provenance  # noqa: E402
from qfab.qasm import QASM  # noqa: E402

CHIP = sys.argv[1] if len(sys.argv) > 1 else "Baihua"
TODAY = time.strftime("%Y-%m-%d")
RUN = time.strftime("%Y-%m-%d_%H%M")   # 同日多轮互不覆盖
EXP = f"e03_bell/{RUN}"
OUT = Path(__file__).resolve().parent.parent / "results" / "e03_bell"
SHOTS = 2048
# Shenglian 单任务宽度墙:>32~40 比特的任务提交后永不返回
# → 每个颜色类再切块,单任务 ≤14 对 = 28 比特,留安全余量
MAX_PAIRS = {"Shenglian": 14}.get(CHIP, 10**9)


def edge_coloring(edges: list[tuple[int, int]]) -> list[list[tuple[int, int]]]:
    """贪心边着色:同色边不共享比特(度≤3 → 颜色数很小)。"""
    classes: list[list[tuple[int, int]]] = []
    for e in edges:
        for cls in classes:
            if all(not (set(e) & set(x)) for x in cls):
                cls.append(e)
                break
        else:
            classes.append([e])
    return classes


def batch_circuit(pairs: list[tuple[int, int]], basis: str) -> tuple[str, list[int]]:
    n = 2 * len(pairs)
    qc = QASM(n)
    tq = []
    for i, (a, b) in enumerate(pairs):
        qc.h(2 * i).cx(2 * i, 2 * i + 1)
        tq += [a, b]
    if basis == "XX":
        for i in range(n):
            qc.h(i)
    qc.barrier().measure_all()
    return qc.qasm(), tq


def main():
    print(f"后端状态: {client.status()}")
    ch = chipmod.load(CHIP, refresh=True)
    av = set(ch.available_qubits())
    edges = [(a, b) for a, b, f in ch.edges
             if a in av and b in av and (f or 0) > 0.5]
    classes = edge_coloring(edges)
    classes = [cls[i:i + MAX_PAIRS] for cls in classes
               for i in range(0, len(cls), MAX_PAIRS)]
    print(f"{len(edges)} 条可测边 → {len(classes)} 个颜色类(含宽度切块)× 2 基 = "
          f"{2 * len(classes)} 个任务")

    specs = []
    for ci, cls in enumerate(classes):
        for basis in ("ZZ", "XX"):
            circ, tq = batch_circuit(cls, basis)
            specs.append(dict(circuit=circ, chip=CHIP, shots=SHOTS,
                              name=f"bell_c{ci}_{basis}", exp=EXP,
                              label=f"class{ci}_{basis}",
                              options={"target_qubits": tq},
                              meta={"class": ci, "basis": basis,
                                    "pairs": [list(p) for p in cls]}))
    recs = client.run_batch(specs, timeout=1200)

    atlas: dict[str, dict] = {}
    for r in recs:
        c = client.counts(r)
        if not c:
            continue
        pairs = r["meta"]["pairs"]
        basis = r["meta"]["basis"]
        n = 2 * len(pairs)
        # 每对比特的 parity = P(00)+P(11) - P(01)-P(10)(在该测量基下)
        tot = analysis.shots_of(c)
        for i, (a, b) in enumerate(pairs):
            even = 0
            for key, v in c.items():
                bits = analysis.bits_of(key, n)
                if bits[2 * i] == bits[2 * i + 1]:
                    even += v
            k = f"{a}-{b}"
            atlas.setdefault(k, {"edge": [a, b],
                                 "vendor_cz": ch.edge_fid(a, b)})
            atlas[k][basis] = round(2 * even / tot - 1, 4)
            lo, hi = analysis.wilson_ci(even, tot)  # even-parity 比例 → parity 尺度
            atlas[k][basis + "_ci"] = [round(2 * lo - 1, 4), round(2 * hi - 1, 4)]

    rows = []
    for k, v in atlas.items():
        if "ZZ" in v and "XX" in v:
            v["asym"] = round(v["ZZ"] - v["XX"], 4)
            v["bell_mean"] = round((v["ZZ"] + v["XX"]) / 2, 4)
        rows.append(v)
    rows.sort(key=lambda r: -(r.get("bell_mean") or -1))

    prov = provenance.stamp(
        qubits=sorted({q for e in edges for q in e}),
        tids=provenance.tids_of(recs), chip=CHIP,
        extra={"selection": "all_usable_edges (vendor CZ>0.5)",
               "n_edges": len(edges), "n_classes": len(classes),
               "shots_per_task": SHOTS,
               "parity_ci": "Wilson 95% on even-parity fraction, mapped to 2p-1"})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"atlas_{CHIP}_{RUN}.json").write_text(
        json.dumps({"provenance": prov, "rows": rows}, ensure_ascii=False, indent=2))
    good = [r for r in rows if (r.get("bell_mean") or 0) > 0.8]
    bad = [r for r in rows if (r.get("bell_mean") or 1) < 0.6]
    asym = sorted(rows, key=lambda r: -abs(r.get("asym") or 0))[:8]
    print(f"\n===== Bell 图谱({len(rows)} 边)=====")
    print(f"  bell_mean>0.8: {len(good)} 条;<0.6: {len(bad)} 条")
    print("  ZZ/XX 不对称 top(CZ 相位误差候选):")
    for r in asym:
        print(f"    {r['edge']}: ZZ={r.get('ZZ')} XX={r.get('XX')} "
              f"asym={r.get('asym')} vendorCZ={r['vendor_cz']}")
    print(f"已写 {OUT / f'atlas_{CHIP}_{RUN}.json'}")


if __name__ == "__main__":
    main()
