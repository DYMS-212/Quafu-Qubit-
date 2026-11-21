"""芯片标定解析与选比特。

chip_info 关键结构(实测):
  size, basis_gates=['h','rx','ry','rz','cz'], priority_qubits, calibration_time
  qubits_info: {'Q0': {index,T1,T2,fidelity,'readout g_fidelity','readout e_fidelity',...}}
  couplers_info: {'C0': {index,fidelity,qubits_index:[a,b]}}

说明:选比特以当日自测为准(available_qubits 过滤 T1=0 等不可用比特);厂家发布的
     标定作为并列参考、不作评判(自测 assignment matrix 由每日体检 e14 产出)。
"""
from __future__ import annotations

from pathlib import Path
from statistics import median


def _readout_fid(qd: dict):
    vals = [v for k, v in qd.items()
            if "readout" in k and "fidelity" in k and isinstance(v, (int, float))]
    return sum(vals) / len(vals) if vals else None


class Chip:
    def __init__(self, info: dict):
        self.info = info
        self.size = info.get("size")
        self.basis_gates = info.get("basis_gates", [])
        self.priority_qubits = info.get("priority_qubits", [])
        self.calibration_time = info.get("calibration_time")

        self.qubits: dict[int, dict] = {}
        for _, qd in info.get("qubits_info", {}).items():
            idx = qd.get("index")
            self.qubits[idx] = {
                "index": idx, "T1": qd.get("T1"), "T2": qd.get("T2"),
                "gate_fidelity": qd.get("fidelity"),
                "readout_fidelity": _readout_fid(qd),
                "raw": qd,
            }

        self.edges: list[tuple[int, int, float]] = []
        for _, cd in info.get("couplers_info", {}).items():
            qi = cd.get("qubits_index") or []
            if len(qi) == 2:
                self.edges.append((qi[0], qi[1], cd.get("fidelity")))

    # ------------------------------------------------------------- 选择
    def available_qubits(self) -> list[int]:
        """T1>0 且单比特门保真>0,否则 transpiler 拒(如 Baihua Q66/Q111)。"""
        return sorted(i for i, q in self.qubits.items()
                      if (q.get("T1") or 0) > 0 and (q.get("gate_fidelity") or 0) > 0)

    def edge_fid(self, a: int, b: int):
        for x, y, f in self.edges:
            if {x, y} == {a, b}:
                return f
        return None

    def neighbors(self, q: int) -> list[int]:
        out = set()
        for x, y, _ in self.edges:
            if x == q:
                out.add(y)
            elif y == q:
                out.add(x)
        return sorted(out)

    def graph(self):
        import networkx as nx
        g = nx.Graph()
        for idx, q in self.qubits.items():
            g.add_node(idx, **{k: v for k, v in q.items() if k != "raw"})
        for a, b, f in self.edges:
            g.add_edge(a, b, fidelity=f)
        return g

    def summary(self) -> dict:
        rf = [q["readout_fidelity"] for q in self.qubits.values() if q["readout_fidelity"]]
        ef = [f for *_, f in self.edges if f]
        return {
            "size": self.size,
            "n_qubits": len(self.qubits),
            "n_available": len(self.available_qubits()),
            "n_edges": len(self.edges),
            "calibration_time": self.calibration_time,
            "readout_median": round(median(rf), 4) if rf else None,
            "cz_median": round(median(ef), 4) if ef else None,
        }


def load(chip: str = "Baihua", refresh: bool = False) -> Chip:
    from . import client
    return Chip(client.get_chip_info(chip, refresh=refresh))


def selfcal_eta(chip: str = "Baihua") -> tuple[dict[int, float], str]:
    """最新自测 η 表:Baihua=e02 全芯片自校准;其他机=e15 分组自校准合并。
    返回 ({qubit: η}, 来源描述);无自测数据时回退 vendor 读出保真并注明。"""
    import json
    from . import analysis
    R = Path(__file__).resolve().parent.parent / "results"
    try:
        if chip == "Baihua":
            f = sorted((R / "e02_selfcal").glob("assignment_*.json"))[-1]
            rows = json.loads(f.read_text())
            return {r["qubit"]: r["eta"] for r in rows}, f.name
        for d in sorted((R / "e15_crossbackend").glob(f"{chip}_*"), reverse=True):
            pairs: dict[tuple, dict[int, dict]] = {}
            for f in sorted(d.glob("selfcal_prep*__*.json")):
                r = json.loads(f.read_text())
                if not r.get("result", {}).get("count"):
                    continue
                p = 1 if r["meta"]["prep"] == "prep1" else 0
                pairs.setdefault(tuple(r["meta"]["qubits"]), {})[p] = r
            eta: dict[int, float] = {}
            for qs, pr in pairs.items():
                if 0 in pr and 1 in pr:
                    asg = analysis.assignment_from_p1(
                        analysis.bit_p1(pr[0]["result"]["count"], len(qs)),
                        analysis.bit_p1(pr[1]["result"]["count"], len(qs)))
                    eta.update({q: asg[i]["eta"] for i, q in enumerate(qs)})
            if eta:
                return eta, d.name
    except Exception:  # noqa: BLE001
        pass
    ch = load(chip)
    return ({q: (v.get("readout_fidelity") or 0) for q, v in ch.qubits.items()},
            "vendor_readout(无自测数据回退)")
