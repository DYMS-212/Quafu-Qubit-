"""e24 · 逻辑 Bell 态:两个逻辑比特之间的纠缠(transversal 逻辑 CZ)。

物理基础(statevector 已验证):[[4,2,2]] 的 transversal S(S⊗4)
在逻辑基下 = diag(1,-1,-1,-1) = **逻辑 CZ × Z̄₁Z̄₂**——一个纠缠性的逻辑双比特门。
流程:|0̄0̄⟩(星形 GHZ)→ H⊗4(|+̄+̄⟩)→ S⊗4 → 逻辑 Bell 类态。

产出态 = 逻辑 cluster 态(Bell 类,单边 H̄ 等价于逻辑 Bell),其三个稳定子
X̄₁Z̄₂ = −Y₀Y₁、Z̄₁X̄₂ = −Y₀Y₂、Ȳ₁Ȳ₂ ∝ Y₁Y₂ **连同检错稳定子 YYYY 全部
落在同一个物理 Y⊗4 读出基**(每比特 S†H 后测)——一个设置测完全部,
F = (1 + Σ_k i_k·c_k)/4 对 Bell 类态是精确值;跑两轮独立统计。
后选 Y⊗4 总奇偶(单比特 X/Z 错均反对易于 YYYY → 可检)。
理想符号 i_k 与后选奇偶由 statevector 现场自校准(绝不手推符号)。

对照臂(同批):物理 Bell,ZZ/XX/YY 三设置。headline 问题:**编码层的纠缠
vs 物理层的纠缠,谁的保真度高?** 纠缠见证:i₁c₁ + i₂c₂ > 1。

运行: ./qrun experiments/e24_logical_bell.py [chip]   # 默认 Baihua;三机全定点
      ./qrun experiments/e24_logical_bell.py selftest  # 零额度自检
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, client, provenance  # noqa: E402
from qfab.qasm import QASM  # noqa: E402

PROJ = Path(__file__).resolve().parent.parent
CHIP = sys.argv[1] if len(sys.argv) > 1 else "Baihua"
RUN = time.strftime("%Y-%m-%d_%H%M")
EXP = f"e24_logical_bell/{RUN}"
OUT = PROJ / "results" / "e24_logical_bell"
SHOTS = 4096
SETTINGS = ("ZZ", "XX", "YY")            # 物理臂三设置
LCORRS = ("c01", "c02", "c12")           # 逻辑臂三关联子(Y⊗4 基,b_i⊕b_j)


# ------------------------------------------------------------------ 线路
def circ_logical() -> str:
    qc = QASM(4)
    qc.h(0).cx(0, 1).cx(0, 2).cx(0, 3)      # |00>_L
    qc.barrier()
    for q in range(4):
        qc.h(q)                              # → |++>_L(含逻辑 SWAP 重标)
    qc.barrier()
    for q in range(4):
        qc.s(q)                              # transversal S = 逻辑 CZ·Z̄₁Z̄₂
    qc.barrier()
    for q in range(4):                       # Y⊗4 读出基(S†H)
        qc.sdg(q)
    for q in range(4):
        qc.h(q)
    qc.barrier().measure_all()
    return qc.qasm()


def circ_physical(setting: str) -> str:
    qc = QASM(2)
    qc.h(0).cx(0, 1)
    qc.barrier()
    if setting == "XX":
        qc.h(0).h(1)
    elif setting == "YY":
        qc.sdg(0).sdg(1).h(0).h(1)
    qc.barrier().measure_all()
    return qc.qasm()


# ------------------------------------------------------------------ 提取
def corr_logical(count: dict, ps_parity: int = 0) -> tuple[dict, float, int]:
    """Y⊗4 基:后选总奇偶 == ps_parity → 三个关联子 (-1)^{b_i⊕b_j}。
    返回 ({c01,c02,c12}, acceptance, n_acc)。"""
    tot = acc = 0
    s = {k: 0 for k in LCORRS}
    for key, v in count.items():
        b = analysis.bits_of(key, 4)
        tot += v
        if (b[0] ^ b[1] ^ b[2] ^ b[3]) == ps_parity:
            acc += v
            s["c01"] += v * (1 - 2 * (b[0] ^ b[1]))
            s["c02"] += v * (1 - 2 * (b[0] ^ b[2]))
            s["c12"] += v * (1 - 2 * (b[1] ^ b[2]))
    return ({k: (x / acc if acc else 0.0) for k, x in s.items()},
            acc / (tot or 1), acc)


def corr_physical(count: dict) -> float:
    tot = s = 0
    for key, v in count.items():
        b = analysis.bits_of(key, 2)
        tot += v
        s += v * (1 - 2 * (b[0] ^ b[1]))
    return s / (tot or 1)


def ideal_signs() -> dict:
    """statevector 现场自校准:后选奇偶 + 全部理想符号(与线路永远同步)。"""
    from tools.verify_circuits_offline import simulate
    cnt = simulate(circ_logical())
    parity = None
    for p in (0, 1):
        _, a, _ = corr_logical(cnt, p)
        if a > 1 - 1e-9:
            parity = p
    assert parity is not None, "码空间 Y⊗4 奇偶不定,线路有误"
    corrs, _, _ = corr_logical(cnt, parity)
    lsigns = {}
    for k, c in corrs.items():
        assert abs(abs(c) - 1) < 1e-9, f"logical/{k} 理想关联子 |{c}| ≠ 1"
        lsigns[k] = round(c)
    psigns = {}
    for st_ in SETTINGS:
        c = corr_physical(simulate(circ_physical(st_)))
        assert abs(abs(c) - 1) < 1e-9, f"physical/{st_} 理想关联子 |{c}| ≠ 1"
        psigns[st_] = round(c)
    return {"logical": lsigns, "physical": psigns, "ps_parity": parity}


def fidelity(corrs: dict, signs: dict) -> float:
    return (1 + sum(signs[k] * corrs[k] for k in signs)) / 4


# ------------------------------------------------------------------ 主流程
def pick_sites():
    """星=度3顶点按瓶颈边(第三好)选;对=自测图谱最优,缺则厂商最优。三机全定点。"""
    from qfab import anchors as anchmod, chip as chipmod
    star = None
    ch = av = None
    if CHIP == "Baihua":
        st_ = anchmod.latest(CHIP).get("star") or {}
        if st_.get("leaves"):
            star = [st_["center"], *st_["leaves"]]
    else:
        ch = chipmod.load(CHIP, refresh=True)
        av = set(ch.available_qubits())
        best = None
        for c in av:
            nb = sorted(((ch.edge_fid(c, x) or 0, x) for x in ch.neighbors(c)
                         if x in av and (ch.edge_fid(c, x) or 0) > 0.5),
                        reverse=True)
            if len(nb) >= 3 and (best is None or nb[2][0] > best[0]):
                best = (nb[2][0], [c] + [x for _, x in nb[:3]])
        if best:
            star = best[1]
    pair = None
    atlas = sorted((PROJ / "results" / "e03_bell").glob(f"atlas_{CHIP}_*.json"))
    if atlas:
        rows = json.loads(atlas[-1].read_text()).get("rows", [])
        best = max((r for r in rows if r.get("bell_mean")),
                   key=lambda r: r["bell_mean"], default=None)
        if best:
            pair = list(best["edge"])
    if not pair and ch is not None:
        edges = [(f, a, b) for a, b, f in ch.edges
                 if a in av and b in av and (f or 0) > 0.5]
        if edges:
            pair = list(max(edges)[1:])
    return star, pair


def selftest():
    signs = ideal_signs()
    for arm in ("logical", "physical"):
        f = fidelity({k: float(v) for k, v in signs[arm].items()}, signs[arm])
        ok = abs(f - 1) < 1e-9
        print(f"  {'✅' if ok else '❌'} {arm}: 理想符号 {signs[arm]} → F={f:.6f}"
              + (f"(后选奇偶={signs['ps_parity']})" if arm == "logical" else ""))
        if not ok:
            sys.exit(1)
    print("selftest 全部通过 ✅")


def main():
    if CHIP == "selftest":
        selftest()
        return
    st = client.status()
    print(f"后端状态: {st}")
    if not isinstance(st.get(CHIP), int):
        print(f"⚠️ {CHIP} 当前不可用({st.get(CHIP)}),中止。")
        return
    signs = ideal_signs()
    star, pair = pick_sites()
    print(f"选址: 逻辑星={star} 物理对={pair}")

    specs = []
    for r_ in (0, 1):                       # 逻辑臂:同一线路两轮独立统计
        spec = dict(circuit=circ_logical(), chip=CHIP, shots=SHOTS,
                    name=f"lbell_L_r{r_}", exp=EXP, label=f"logical_r{r_}",
                    meta={"arm": "logical", "rep": r_})
        if star:
            spec["options"] = {"target_qubits": star}
        specs.append(spec)
    for s in SETTINGS:
        spec = dict(circuit=circ_physical(s), chip=CHIP, shots=SHOTS,
                    name=f"lbell_P_{s}", exp=EXP, label=f"physical_{s}",
                    meta={"arm": "physical", "setting": s})
        if pair:
            spec["options"] = {"target_qubits": pair}
        specs.append(spec)
    print(f"共 {len(specs)} 条线路 × {SHOTS} shots")
    recs = client.run_batch(specs, timeout=1200)

    rows = []
    summary_arms = {}
    # 逻辑臂:合并两轮 counts
    merged: dict = {}
    for r in recs:
        if r["meta"]["arm"] == "logical" and client.counts(r):
            for k, v in client.counts(r).items():
                merged[k] = merged.get(k, 0) + v
    if merged:
        corrs, a, n_acc = corr_logical(merged, signs["ps_parity"])
        f = fidelity(corrs, signs["logical"])
        wit = sum(signs["logical"][k] * corrs[k] for k in ("c01", "c02"))
        summary_arms["logical"] = {
            "fidelity": round(f, 5), "witness": round(wit, 5),
            "entangled_witness": bool(wit > 1), "acceptance": round(a, 4),
            "n_accepted": n_acc,
            "corrs": {k: round(v, 5) for k, v in corrs.items()}}
        rows += [{"arm": "logical", "corr_name": k, "corr": round(v, 5),
                  "ideal_sign": signs["logical"][k], "acceptance": round(a, 4)}
                 for k, v in corrs.items()]
        print(f"  logical:  F={f:.4f} 见证={wit:.3f}"
              f"({'>1 逻辑纠缠 ✅' if wit > 1 else '≤1 未见证'}) acc={a:.3f}")
    # 物理臂
    pcorrs = {}
    for r in recs:
        if r["meta"]["arm"] != "physical" or not client.counts(r):
            continue
        s = r["meta"]["setting"]
        pcorrs[s] = corr_physical(client.counts(r))
        rows.append({"arm": "physical", "corr_name": s,
                     "corr": round(pcorrs[s], 5), "ideal_sign": signs["physical"][s],
                     "acceptance": 1.0})
    if len(pcorrs) == 3:
        f = fidelity(pcorrs, signs["physical"])
        wit = sum(signs["physical"][s] * pcorrs[s] for s in ("ZZ", "XX"))
        summary_arms["physical"] = {
            "fidelity": round(f, 5), "witness": round(wit, 5),
            "entangled_witness": bool(wit > 1),
            "corrs": {k: round(v, 5) for k, v in pcorrs.items()}}
        print(f"  physical: F={f:.4f} 见证={wit:.3f}"
              f"({'>1 纠缠 ✅' if wit > 1 else '≤1 未见证'})")

    OUT.mkdir(parents=True, exist_ok=True)
    summary = {"provenance": provenance.stamp(
                   anchor=None, chip=CHIP, tids=provenance.tids_of(recs),
                   qubits={"logical_star": star, "physical_pair": pair},
                   extra={"gate": "transversal S = logical CZ·Z1Z2 (statevector 已验证)",
                          "ideal_signs": signs}),
               "chip": CHIP, "arms": summary_arms, "rows": rows}
    (OUT / f"summary_{CHIP}_{RUN}.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"已写 {OUT / f'summary_{CHIP}_{RUN}.json'}")


if __name__ == "__main__":
    main()
