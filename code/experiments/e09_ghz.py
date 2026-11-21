"""e09 · GHZ 纠缠规模:线性链 vs 拓扑感知树形(depth~√N),扫到断点。

指标:Z 基布居 P(全0)+P(全1);X 基相干 ⟨X^N⟩(读出平均伪造不了的真量子量)。
线性链:锚点链前缀(N≤13,对照用);树形:最优比特 BFS 生成树(父先子后 CX),
N 扫到芯片/宽度墙允许的最大——纠缠规模在哪个 N 断、断在哪台机,是三机能力
地图最直观的一行。

运行: ./qrun experiments/e09_ghz.py [chip]   # 默认 Baihua
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, anchors as anchmod, chip as chipmod, client, provenance  # noqa: E402
from qfab.qasm import QASM  # noqa: E402

CHIP = sys.argv[1] if len(sys.argv) > 1 else "Baihua"
TODAY = time.strftime("%Y-%m-%d")
RUN = time.strftime("%Y-%m-%d_%H%M")   # 同日多轮互不覆盖
TAG = RUN if CHIP == "Baihua" else f"{CHIP}_{RUN}"
EXP = f"e09_ghz/{TAG}"
PROJ = Path(__file__).resolve().parent.parent
OUT = PROJ / "results" / "e09_ghz"
SHOTS = 4096
# 树形扫到断点;线性只做 ≤13 的布线对照。Shenglian 受单任务宽度墙(32-40)限制 N≤25。
NS_TREE = {"Shenglian": [3, 5, 7, 9, 11, 13, 17, 21, 25]}.get(
    CHIP, [3, 5, 7, 9, 11, 13, 17, 21, 25, 31])
NS_LIN = [3, 5, 7, 9, 11, 13]


def load_eta() -> tuple[dict[int, float], str]:
    return chipmod.selfcal_eta(CHIP)


def linear_chain() -> list[int]:
    """线性布线对照的物理链:锚点最长链(选址与当日体检同源)。"""
    ch_ = anchmod.latest(CHIP).get("chains") or {}
    longest = max(ch_.values(), key=len, default=[])
    return list(longest)


def bfs_tree(ch, eta, n: int):
    """从最优比特 BFS 生成 n 节点树。返回 (节点序, 边=(父idx,子idx) 列表)。"""
    av = set(ch.available_qubits()) & {q for q, e in eta.items() if e > 0.4}
    root = max(av, key=lambda q: (len([x for x in ch.neighbors(q) if x in av]), eta[q]))
    order, edges, seen = [root], [], {root}
    frontier = [root]
    while frontier and len(order) < n:
        nxt = []
        for p in frontier:
            for x in sorted((x for x in ch.neighbors(p)
                             if x in av and x not in seen
                             and (ch.edge_fid(p, x) or 0) > 0.5),
                            key=lambda x: -eta[x]):
                if len(order) >= n:
                    break
                seen.add(x)
                edges.append((order.index(p), len(order)))
                order.append(x)
                nxt.append(x)
        frontier = nxt
    return order, edges


def ghz_circuit(n: int, edges: list[tuple[int, int]], basis: str) -> str:
    qc = QASM(n)
    qc.h(0)
    for p, c in edges:            # 父先子后(BFS 序保证)
        qc.cx(p, c)
    if basis == "X":
        for i in range(n):
            qc.h(i)
    qc.barrier().measure_all()
    return qc.qasm()


def main():
    print(f"后端状态: {client.status()}")
    ch = chipmod.load(CHIP, refresh=True)
    eta, eta_src = load_eta()
    lin = linear_chain()
    specs = []
    for n in sorted(set(NS_TREE) | set(NS_LIN)):
        lin_edges = [(i, i + 1) for i in range(n - 1)]
        tree_nodes, tree_edges = bfs_tree(ch, eta, n)
        layouts = {}
        if n in NS_LIN and len(lin) >= n:
            layouts["linear"] = (lin[:n], lin_edges)
        if n in NS_TREE:
            layouts["tree"] = (tree_nodes, tree_edges)
        for name, (tq, edges) in layouts.items():
            if len(tq) < n:
                print(f"跳过 {name} N={n}(可用优质子图不足)")
                continue
            depth = max((0, *[_depth(edges, i) for i in range(1, n)]))
            for basis in ("Z", "X"):
                specs.append(dict(circuit=ghz_circuit(n, edges, basis),
                                  chip=CHIP, shots=SHOTS,
                                  name=f"ghz_{name}_N{n}_{basis}", exp=EXP,
                                  label=f"{name}_N{n}_{basis}",
                                  options={"target_qubits": tq},
                                  meta={"layout": name, "N": n, "basis": basis,
                                        "qubits": tq, "depth": depth}))
    print(f"共 {len(specs)} 条线路")
    recs = client.run_batch(specs, timeout=1200)

    rows = []
    for r in recs:
        c = client.counts(r)
        if not c:
            continue
        m = r["meta"]
        n = m["N"]
        tot = analysis.shots_of(c)
        if m["basis"] == "Z":
            pop_ct = sum(v for k, v in c.items()
                         if set(k.zfill(n)) <= {"0"} or set(k.zfill(n)) <= {"1"})
            lo, hi = analysis.wilson_ci(pop_ct, tot)
            rows.append({**m, "chip": CHIP, "metric": "pop00_11", "value": round(pop_ct / tot, 4),
                         "ci": [round(lo, 4), round(hi, 4)]})
        else:
            even = sum(v for k, v in c.items() if k.count("1") % 2 == 0)
            lo, hi = analysis.wilson_ci(even, tot)  # even 比例 → ⟨X^N⟩ 尺度
            rows.append({**m, "chip": CHIP, "metric": "X_parity",
                         "value": round(2 * even / tot - 1, 4),
                         "ci": [round(2 * lo - 1, 4), round(2 * hi - 1, 4)]})
    prov = provenance.stamp(
        qubits=None, tids=provenance.tids_of(recs), chip=CHIP,
        extra={"selection": "linear=锚点链前缀; tree=最优比特 BFS(η>0.4, edge CZ>0.5)",
               "linear_chain": lin, "eta_source": eta_src,
               "Ns_tree": NS_TREE, "Ns_linear": NS_LIN,
               "shots_per_task": SHOTS,
               "note": "每行物理比特见 row['qubits']",
               "ci": "Wilson 95%: pop 为比例;X_parity 映射到 2p-1"})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"summary_{TAG}.json").write_text(
        json.dumps({"provenance": prov, "rows": rows}, ensure_ascii=False, indent=2))
    print("\n===== GHZ 规模(pop=布居, X=相干)=====")
    for n in sorted(set(NS_TREE) | set(NS_LIN)):
        line = f"  N={n:2d}: "
        for lay in ("linear", "tree"):
            pop = [r["value"] for r in rows if r["layout"] == lay
                   and r["N"] == n and r["metric"] == "pop00_11"]
            xp = [r["value"] for r in rows if r["layout"] == lay
                  and r["N"] == n and r["metric"] == "X_parity"]
            if pop and xp:
                line += f"{lay}: pop={pop[0]:.3f} ⟨X^N⟩={xp[0]:+.3f}   "
        print(line)
    print(f"已写 {OUT / f'summary_{TAG}.json'}")


def _depth(edges: list[tuple[int, int]], node: int) -> int:
    d = 0
    cur = node
    parents = {c: p for p, c in edges}
    while cur in parents:
        cur = parents[cur]
        d += 1
    return d


if __name__ == "__main__":
    main()
