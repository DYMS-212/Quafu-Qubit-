"""每机锚点管理:首次自动生成 → 冻结到 results/anchors/<chip>.json(跨日可比的命根子)。

锚点内容:rep 链(d=3..13 嵌套前缀)、syndrome 三元组、Bell 边、负载边、物理基线比特、最优星形。
生成依据:当机自测 η(若有 e14/e15 selfcal)> vendor 读出;边一律要求 vendor CZ>0.5。
重新生成:delete 配置文件或 load(chip, force=True)——注意这会破坏该机的跨日可比性。
"""
from __future__ import annotations

import json
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
ANCHOR_DIR = PROJ / "results" / "anchors"


def _usable_graph(ch, min_edge: float = 0.5):
    av = set(ch.available_qubits())
    adj: dict[int, list[int]] = {q: [] for q in av}
    for a, b, f in ch.edges:
        if a in av and b in av and (f or 0) > min_edge:
            adj[a].append(b)
            adj[b].append(a)
    return av, adj


def _score_map(ch, eta: dict[int, float] | None, mode: str = "readout"):
    """比特评分。mode:
      readout   — 读出 η(有自测用 η,否则 vendor 读出);旧默认,保持跨日可比。
      gate      — 单比特门保真(vendor)。
      composite — 读出 × 门 的几何平均(兼顾读出与门,供编码/深门电路选址)。
    """
    ro = dict(eta) if eta else {q: (d.get("readout_fidelity") or 0)
                                for q, d in ch.qubits.items()}
    if mode == "readout":
        return ro
    gf = {q: (d.get("gate_fidelity") or 0) for q, d in ch.qubits.items()}
    if mode == "gate":
        return gf
    return {q: (max(ro.get(q, 0), 0) * max(gf.get(q, 0), 0)) ** 0.5
            for q in set(ro) | set(gf)}


def _best_chain(adj, score, n: int, edge_score=None):
    """长度 n 的最优路径。node 分 = min(score);若给 edge_score(a,b)→CZ 质量,
    则路径分 = min(node) × min(沿途 CZ),让链选址也吃门质量。"""
    best = None
    for start in adj:
        stack = [[start]]
        while stack:
            p = stack.pop()
            if len(p) == n:
                s = min(score.get(q, 0) for q in p)
                if edge_score is not None and len(p) > 1:
                    s *= min((edge_score(p[i], p[i + 1]) or 0) for i in range(len(p) - 1))
                if best is None or s > best[0]:
                    best = (s, p)
                continue
            for x in adj[p[-1]]:
                if x not in p:
                    stack.append(p + [x])
    return best[1] if best else None


def _long_chain(adj, score, n_max: int, edge_score=None, beam: int = 400):
    """束搜索最长优质路径(≤n_max 比特)。全枚举 DFS 在 25 比特长链上会爆炸,
    束搜索按 分数=min(节点)×min(边) 保留 top-beam 局部路径,有界且够优。
    返回能达到的最长长度下的最优路径(None=图为空)。"""
    frontier = sorted((((score.get(q, 0), 1.0), [q]) for q in adj),
                      key=lambda t: -(t[0][0] * t[0][1]))[:beam]
    best_at: dict[int, list[int]] = {1: frontier[0][1]} if frontier else {}
    length = 1
    while frontier and length < n_max:
        nxt = []
        for (ns, es), p in frontier:
            for x in adj[p[-1]]:
                if x in p:
                    continue
                ns2 = min(ns, score.get(x, 0))
                es2 = es if edge_score is None else min(
                    es, (edge_score(p[-1], x) or 0))
                nxt.append(((ns2, es2), p + [x]))
        if not nxt:
            break
        nxt.sort(key=lambda t: -(t[0][0] * t[0][1]))
        frontier = nxt[:beam]
        length += 1
        best_at[length] = frontier[0][1]
    for ln in range(n_max, 0, -1):
        if best_at.get(ln):
            return best_at[ln]
    return None


def _generate(ch, eta: dict[int, float] | None, mode: str = "readout") -> dict:
    av, adj = _usable_graph(ch)
    score = _score_map(ch, eta, mode)
    es = None if mode == "readout" else ch.edge_fid  # 门感知模式:链选址吃 CZ 质量
    anchors: dict = {"chip_calibration_time": ch.calibration_time,
                     "generated_from": "selfcal_eta" if eta else "vendor_readout",
                     "score_mode": mode}
    # 链(嵌套:d13 ⊃ d11 ⊃ … ⊃ d3,前缀共用,便于同区域对比;
    # 束搜索找 ≤25 比特最长优质路径,d3 前缀更敏感 → 取 d3 段更好的方向)
    anchors["chains"] = {}
    c_long = _long_chain(adj, score, 25, es)
    if c_long and len(c_long) >= 5:
        def _seg_score(p):
            s = min(score.get(q, 0) for q in p)
            if es is not None:
                s *= min((es(p[i], p[i + 1]) or 0) for i in range(len(p) - 1))
            return s
        if _seg_score(c_long[::-1][:5]) > _seg_score(c_long[:5]):
            c_long = c_long[::-1]
        anchors["chains"] = {str(d): c_long[:2 * d - 1]
                             for d in (3, 5, 7, 9, 11, 13)
                             if 2 * d - 1 <= len(c_long)}
    # syndrome 三元组:最优 3 比特路径 [D0, A, D1]
    t = _best_chain(adj, score, 3, es)
    if t:
        anchors["syn_triple"] = {"data": [t[0], t[2]], "anc": t[1]}
    # Bell 边:两条按 min(score) 最好的边(不相交)
    edges = sorted(((min(score.get(a, 0), score.get(b, 0)), a, b)
                    for a in adj for b in adj[a] if a < b), reverse=True)
    bell, used = [], set()
    for _, a, b in edges:
        if not ({a, b} & used):
            bell.append([a, b])
            used |= {a, b}
        if len(bell) == 2:
            break
    anchors["bell_edges"] = bell
    # 负载边:离主链最远的可用边(粗略:不与链/bell 重叠即可)
    chain_set = set(anchors["chains"].get("7", anchors["chains"].get("3", [])))
    load = next(([a, b] for _, a, b in edges
                 if not ({a, b} & chain_set) and not ({a, b} & used)), None)
    anchors["load_edge"] = load
    # 物理基线:score 最高的 2 个比特
    anchors["phys_baseline"] = sorted(av, key=lambda q: -score.get(q, 0))[:2]
    # 最优星形(度3)
    best_star = None
    for c in av:
        nb = [x for x in adj[c]]
        if len(nb) >= 3:
            leaves = sorted(nb, key=lambda x: -score.get(x, 0))[:3]
            s = min(score.get(q, 0) for q in [c, *leaves])
            if best_star is None or s > best_star[0]:
                best_star = (s, c, leaves)
    if best_star:
        anchors["star"] = {"center": best_star[1], "leaves": best_star[2]}
    return anchors


def assignment_from_probe_dir(day: Path):
    """从一个体检日目录选最新的**有效** selfcal_prep0/prep1 对,重算读出 assignment。

    队列超时未补收的任务 result=null → 跳过;文件名含 tid(时间戳)倒序取最新有效者,
    当天多轮体检时不因最后一轮尚未补收就整体放弃(选址在队列堵塞下仍稳)。
    返回 (qubits, asg);无任何有效对时返回 (None, None)。
    """
    from . import analysis

    def _latest_valid(tag: str):
        for f in sorted(day.glob(f"{tag}__*.json"), reverse=True):
            rec = json.loads(f.read_text())
            if (rec.get("result") or {}).get("count"):
                return rec
        return None

    rec0, rec1 = _latest_valid("selfcal_prep0"), _latest_valid("selfcal_prep1")
    if not (rec0 and rec1):
        return None, None
    qubits = rec0["meta"]["qubits"]
    n = len(qubits)
    asg = analysis.assignment_from_p1(
        analysis.bit_p1(rec0["result"]["count"], n),
        analysis.bit_p1(rec1["result"]["count"], n))
    return qubits, asg


def latest_probe_assignment(probe_root: Path):
    """取最近一天有**有效**自校准的体检 assignment(离线分析用)。

    今天体检若全超时未补收,自动回退到最近一个有有效数据的日子(打印实际用的日期,透明可溯源)。
    返回 (qubits, asg, 体检日目录名);全无有效时 (None, None, None)。
    """
    if not probe_root.exists():
        return None, None, None
    days = sorted((d for d in probe_root.iterdir() if d.is_dir()), reverse=True)
    for day in days:
        qubits, asg = assignment_from_probe_dir(day)
        if qubits:
            return qubits, asg, day.name
    return None, None, None


def daily_best(chip_name: str, eta: dict[int, float], date: str,
               mode: str = "readout") -> dict:
    """每日自适应选址:用当天自测 η 重算最优锚点,落盘留痕(不动冻结锚点)。
    双层机制:e14 体检用冻结锚点(漂移标尺);headline 实验用本函数的当日最优。
    mode: readout(默认,保跨日可比)/ gate / composite(门感知,e21 验证用)。"""
    from . import chip as chipmod
    ch = chipmod.load(chip_name, refresh=True)
    anchors = _generate(ch, eta, mode)
    anchors["date"] = date
    d = ANCHOR_DIR / "daily"
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{chip_name}_{date}.json").write_text(
        json.dumps(anchors, ensure_ascii=False, indent=2))
    return anchors


def latest(chip_name: str) -> dict:
    """headline 实验选址入口:最新的当日最优锚点;没有则退回冻结锚点。"""
    d = ANCHOR_DIR / "daily"
    if d.exists():
        files = sorted(d.glob(f"{chip_name}_*.json"))
        if files:
            a = json.loads(files[-1].read_text())
            print(f"选址: 当日最优锚点 {files[-1].name}")
            return a
    print("选址: 冻结锚点(当日最优不存在,先跑 e14+e02)")
    return load(chip_name)


def load(chip_name: str, eta: dict[int, float] | None = None,
         force: bool = False) -> dict:
    """取该机冻结锚点;不存在(或 force)则由当前标定生成并冻结。"""
    ANCHOR_DIR.mkdir(parents=True, exist_ok=True)
    f = ANCHOR_DIR / f"{chip_name}.json"
    if f.exists() and not force:
        return json.loads(f.read_text())
    from . import chip as chipmod
    ch = chipmod.load(chip_name, refresh=True)
    anchors = _generate(ch, eta)
    f.write_text(json.dumps(anchors, ensure_ascii=False, indent=2))
    print(f"⚠️ {chip_name} 锚点已生成并冻结 → {f}")
    return anchors
