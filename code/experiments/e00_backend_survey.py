"""e00 · 后端普查 + 超参兼容性审计(纯读标定,零额度)。

对每个后端:size/比特数/可用数/度分布/标定边/读出与CZ中位/basis gates;
并审计本套件的 Baihua 锚点超参(链/锚/边/基线比特)在该机是否成立。

运行: ./qrun experiments/e00_backend_survey.py
"""
from __future__ import annotations

import json
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import chip as chipmod, client, provenance  # noqa: E402

TODAY = time.strftime("%Y-%m-%d")
OUT = Path(__file__).resolve().parent.parent / "results" / "e00_survey"

# 本套件当前写死的 Baihua 锚点(audit 对象)
ANCHORS = {
    "rep_chain_d3": [8, 7, 6, 5, 4],
    "rep_chain_d7": [8, 7, 6, 5, 4, 3, 2, 15, 14, 13, 26, 27, 28],
    "syn_triple": [20, 22, 21],
    "bell_edges": [(67, 68), (73, 74)],
    "load_edge": [(73, 74)],
    "phys_baseline": [8, 150],
}


def audit_chip(ch: chipmod.Chip) -> dict:
    av = set(ch.available_qubits())

    def chain_ok(chain):
        if not all(q in av for q in chain):
            return "缺比特"
        for a, b in zip(chain, chain[1:]):
            if (ch.edge_fid(a, b) or 0) <= 0.5:
                return f"边 {a}-{b} 不可用"
        return "OK"

    def edges_ok(edges):
        bad = [f"{a}-{b}" for a, b in edges
               if a not in av or b not in av or (ch.edge_fid(a, b) or 0) <= 0.5]
        return "OK" if not bad else f"不可用: {bad}"

    return {
        "rep_chain_d3": chain_ok(ANCHORS["rep_chain_d3"]),
        "rep_chain_d7": chain_ok(ANCHORS["rep_chain_d7"]),
        "syn_triple(20-21,22-21)": edges_ok([(20, 21), (22, 21)]),
        "bell_edges": edges_ok(ANCHORS["bell_edges"]),
        "load_edge": edges_ok(ANCHORS["load_edge"]),
        "phys_baseline": "OK" if all(q in av for q in ANCHORS["phys_baseline"])
                          else "缺比特",
    }


def main():
    st = client.status()
    print(f"后端状态: {st}\n")
    report = {"date": TODAY, "status": st, "chips": {}}
    # 溯源戳:零额度纯读普查,无 shots/无比例指标(故无 CI);记录审计的是哪套锚点、扫了哪些后端。
    report["provenance"] = provenance.stamp(
        anchor=None, qubits=None, tids=None, chip="(survey)",
        extra={"anchor": "frozen_hardcoded",
               "kind": "backend_survey_zero_quota",
               "audited_anchor_set": "Baihua_hardcoded",
               "backends": list(st.keys())})
    for name, q in st.items():
        entry: dict = {"queue": q, "usable_now": isinstance(q, int)}
        try:
            ch = chipmod.load(name, refresh=True)
            av = ch.available_qubits()
            deg = Counter(len([x for x in ch.neighbors(qq) if x in set(av)])
                          for qq in av)
            n_edges_ok = sum(1 for a, b, f in ch.edges
                             if a in set(av) and b in set(av) and (f or 0) > 0.5)
            entry.update({
                **ch.summary(),
                "available_n": len(av),
                "degree_distribution": dict(sorted(deg.items())),
                "edges_usable": n_edges_ok,
                "max_chain_needed_d7": 13,
                "anchor_audit": audit_chip(ch),
            })
        except Exception as e:  # noqa: BLE001
            entry["chip_info_error"] = str(e)[:150]
        report["chips"][name] = entry
        print(f"===== {name} =====")
        print(json.dumps(entry, ensure_ascii=False, indent=2))
        print()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"survey_{TODAY}.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2))
    print(f"已写 {OUT / f'survey_{TODAY}.json'}")


if __name__ == "__main__":
    main()
