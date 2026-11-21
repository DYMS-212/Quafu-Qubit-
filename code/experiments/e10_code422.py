"""e10 · [[4,2,2]] 检错码:逻辑比特 vs 物理比特正面对决(headline)。

码:stabilizer S_X=XXXX, S_Z=ZZZZ;逻辑算符(码比特 1-4 ↔ q0-q3):
    X̄₁=X₁X₂  Z̄₁=Z₁Z₃  |  X̄₂=X₁X₃  Z̄₂=Z₁Z₂   (对易/反对易已验证)
  |00⟩_L = (|0000⟩+|1111⟩)/√2 = 4-GHZ → heavy-hex 度3星形(中心+3叶)完美放置
  |++⟩_L = H⊗4|00⟩_L(已验证仍在码空间)

实验矩阵(每个星形 8 条线路):
  Z 基:4 逻辑态 |mn⟩_L(GHZ + 转横 X̄),测 Z⊗4。
       后选 S_Z:总奇偶=0(检出任意单 X 错);读 Z̄₁=b0⊕b2, Z̄₂=b0⊕b1。
  X 基:4 逻辑态 |±±⟩_L(GHZ + H⊗4 + 转横 Z̄,barrier 防 H 对消),测 H⊗4 后 Z。
       后选 S_X;读 X̄₁=b0⊕b1, X̄₂=b0⊕b2(检出任意单 Z 错)。
对照:全芯片自测 η 最好的 2 个物理比特,|0>/|1>/|+>/|−> 裸存读——最严苛基线。
选址:当日体检自测 η(153 比特)对全部度3星打分(min-η 保守分),取 top-N 不相交。

**融合制备臂(验证制备)**:星形 GHZ 制备存在单故障→不可检
weight-2 X(≡逻辑错)的位置(如 cx(0,3) 前后的中心 X 错 → X₀X₃ ≡ X̄₁X̄₂),
这是 e10 距最强物理比特差距的候选来源。度 3 拓扑放不下"星形+验证比特"
(需与两叶同时相邻=四元环,heavy-hex 无),但**线形 5 比特可以**:
  path q0-q1-q2-q3-q4,码比特 (c1..c4)=(q0,q1,q3,q4),q2=验证比特 v。
  Bell(q0,q1) ⊗ Bell(q4,q3) → cx(q1,v), cx(q3,v) → v = b1⊕b3。
  末端后选 v=0 即把 Bell⊗Bell 投影成 4-GHZ(延迟测量原理;理想接受率 1/2)。
**容错账本(statevector 验证器逐故障核对,诚实)**:码比特/Bell 制备/CX 的全部
单 X 故障要么翻 v(验证剔除)、要么单 X(S_Z 剔除)——星形的"中心跨 CX 窗"
类致命故障被消除;**残余暴露=v 自身**:v 从首个 CX 起持有分支比特,其间任何
X 翻转(退相干或 1→0 误读,~ε₁(v)/2)→ 选错融合分支 ≡ X̄₁,不可检。
即危险窗从星形的 2 类压缩到 1 类,且该类集中在单个比特 v 上 → 选址时 v 的
读出质量加权优先。代价=接受率减半。与星形臂、物理基线同批交错,
直接实测"验证制备能否翻盘"。

运行: ./qrun experiments/e10_code422.py [n_stars]   # 默认 6 星 + 3 融合路径
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, anchors, chip as chipmod, client, provenance  # noqa: E402
from qfab.qasm import QASM  # noqa: E402

CHIP = "Baihua"
TODAY = time.strftime("%Y-%m-%d")
RUN = time.strftime("%Y-%m-%d_%H%M")   # 同日多轮互不覆盖
EXP = f"e10_code422/{RUN}"
PROJ = Path(__file__).resolve().parent.parent
OUT = PROJ / "results" / "e10_code422"
SHOTS = 4096
N_STARS = int(sys.argv[1]) if len(sys.argv) > 1 else 6
N_PATHS = 3          # 融合制备(FT 验证)5-path 条数


# ------------------------------------------------------------------ 选址
def load_today_eta() -> dict[int, float]:
    """从当日体检的自校准原始数据重算全芯片 η(跳过队列超时未补收的记录)。"""
    d = PROJ / "results" / "e14_probe" / TODAY
    qubits, asg = anchors.assignment_from_probe_dir(d)
    if not qubits:
        raise RuntimeError("缺当日有效自校准数据(体检超时未补收?先 ./qrun tools/recollect.py)")
    return {q: asg[i]["eta"] for i, q in enumerate(qubits)}


def pick_stars(ch: chipmod.Chip, eta: dict[int, float], n_stars: int):
    av = set(ch.available_qubits()) & set(eta)
    stars = []
    for c in av:
        nb = [x for x in ch.neighbors(c) if x in av
              and (ch.edge_fid(c, x) or 0) > 0.5]   # 未标定边(vendor=0)不可用
        if len(nb) >= 3:
            leaves = sorted(nb, key=lambda x: -eta[x])[:3]
            group = [c, *leaves]
            stars.append({
                "center": c, "leaves": leaves,
                "score_min": min(eta[q] for q in group),
                "score_mean": sum(eta[q] for q in group) / 4,
                "cz": [ch.edge_fid(c, l) for l in leaves],
            })
    stars.sort(key=lambda s: (-s["score_min"], -s["score_mean"]))
    picked, used = [], set()
    for s in stars:
        g = {s["center"], *s["leaves"]}
        if g & used:
            continue
        picked.append(s)
        used |= g
        if len(picked) == n_stars:
            break
    return picked


def pick_paths(ch: chipmod.Chip, eta: dict[int, float], n_paths: int):
    """融合制备用 5-path:q0-q1-v-q3-q4,v(中点)=验证比特,须与 q1/q3 相邻。
    打分含 v(它也被测量);同一 5 比特组合只留最优排列;取 top-N 不相交。"""
    av = set(ch.available_qubits()) & set(eta)

    def ok(a, b):
        return (ch.edge_fid(a, b) or 0) > 0.5

    best: dict[frozenset, dict] = {}
    for v in av:
        mids = [x for x in ch.neighbors(v) if x in av and ok(v, x)]
        for a in mids:
            for b in mids:
                if a == b:
                    continue
                for q0 in (x for x in ch.neighbors(a)
                           if x in av and x not in (v, b) and ok(a, x)):
                    for q4 in (x for x in ch.neighbors(b)
                               if x in av and x not in (v, a, q0) and ok(b, x)):
                        path = [q0, a, v, b, q4]
                        key = frozenset(path)
                        # v 持有分支比特→其 X 翻转/误读=不可检逻辑错,η_v 单列加权
                        cand = {"path": path, "vq": v, "eta_v": eta[v],
                                "score_min": min(eta[q] for q in path),
                                "score_mean": sum(eta[q] for q in path) / 5}
                        cur = best.get(key)
                        if not cur or (cand["score_min"], cand["eta_v"]) > \
                                (cur["score_min"], cur["eta_v"]):
                            best[key] = cand
    paths = sorted(best.values(),
                   key=lambda p: (-p["score_min"], -p["eta_v"], -p["score_mean"]))
    picked, used = [], set()
    for p in paths:
        g = set(p["path"])
        if g & used:
            continue
        picked.append(p)
        used |= g
        if len(picked) == n_paths:
            break
    return picked


# ------------------------------------------------------------------ 线路
def circ_422_z(m: int, n: int) -> str:
    qc = QASM(4)
    qc.h(0).cx(0, 1).cx(0, 2).cx(0, 3)     # GHZ = |00>_L
    qc.barrier()
    if m:                                   # X̄₁ = X q0,q1
        qc.x(0).x(1)
    if n:                                   # X̄₂ = X q0,q2
        qc.x(0).x(2)
    qc.barrier().measure_all()
    return qc.qasm()


def circ_422_x(m: int, n: int) -> str:
    qc = QASM(4)
    qc.h(0).cx(0, 1).cx(0, 2).cx(0, 3)
    for q in range(4):
        qc.h(q)                             # → |++>_L
    qc.barrier()                            # 防止与测量 H 对消
    if m:                                   # Z̄₁ = Z q0,q2 (置 X̄₁=−)
        qc.z(0).z(2)
    if n:                                   # Z̄₂ = Z q0,q1 (置 X̄₂=−)
        qc.z(0).z(1)
    qc.barrier()
    for q in range(4):
        qc.h(q)                             # X 基读出
    qc.measure_all()
    return qc.qasm()


def circ_fusion_z(m: int, n: int) -> str:
    """融合制备 |mn⟩_L,Z 基读出。逻辑位 0..4 = path 顺序,v=逻辑位 2。
    码比特 (c1,c2,c3,c4)=(q0,q1,q3,q4):X̄₁=X(q0,q1),X̄₂=X(q0,q3)。"""
    qc = QASM(5)
    qc.h(0).cx(0, 1)                        # Bell(q0,q1)
    qc.h(4).cx(4, 3)                        # Bell(q4,q3)
    qc.barrier()
    qc.cx(1, 2).cx(3, 2)                    # v = b1⊕b3(末端后选 v=0 → GHZ₄)
    qc.barrier()
    if m:
        qc.x(0).x(1)
    if n:
        qc.x(0).x(3)
    qc.barrier().measure_all()
    return qc.qasm()


def circ_fusion_x(m: int, n: int) -> str:
    """融合制备 + H⊗4 → |±±⟩_L,X 基读出。Z̄₁=Z(q0,q3),Z̄₂=Z(q0,q1)。"""
    qc = QASM(5)
    qc.h(0).cx(0, 1)
    qc.h(4).cx(4, 3)
    qc.barrier()
    qc.cx(1, 2).cx(3, 2)
    qc.barrier()
    for q in (0, 1, 3, 4):
        qc.h(q)                             # → |++⟩_L
    qc.barrier()                            # 防止与测量 H 对消
    if m:
        qc.z(0).z(3)
    if n:
        qc.z(0).z(1)
    qc.barrier()
    for q in (0, 1, 3, 4):
        qc.h(q)                             # X 基读出
    qc.measure_all()
    return qc.qasm()


def circ_phys(basis: str, bit: int) -> str:
    qc = QASM(1)
    if basis == "Z":
        if bit:
            qc.x(0)
    else:
        qc.h(0)
        if bit:
            qc.z(0)
        qc.barrier()
        qc.h(0)
    qc.barrier().measure_all()
    return qc.qasm()


# ------------------------------------------------------------------ 分析
def eval_422(count: dict, basis: str, m: int, n: int) -> dict:
    """后选 stabilizer 奇偶,读逻辑算符。Z基: Z̄₁=b0⊕b2,Z̄₂=b0⊕b1;
    X基: X̄₁=b0⊕b1,X̄₂=b0⊕b2。"""
    tot = acc = e1 = e2 = ej = raw_e1 = raw_e2 = 0
    for key, v in count.items():
        b = analysis.bits_of(key, 4)
        tot += v
        if basis == "Z":
            l1, l2 = b[0] ^ b[2], b[0] ^ b[1]
        else:
            l1, l2 = b[0] ^ b[1], b[0] ^ b[2]
        raw_e1 += v * (l1 != m)
        raw_e2 += v * (l2 != n)
        if (b[0] ^ b[1] ^ b[2] ^ b[3]) == 0:      # stabilizer 后选
            acc += v
            e1 += v * (l1 != m)
            e2 += v * (l2 != n)
            ej += v * ((l1 != m) or (l2 != n))
    a = acc or 1
    return {"shots": tot, "acceptance": acc / (tot or 1),
            "acceptance_ci": analysis.wilson_ci(acc, tot or 1),
            "ps_err_L1": e1 / a, "ps_err_L2": e2 / a, "ps_err_joint": ej / a,
            "ps_err_perq": (e1 + e2) / (2 * a),
            "ps_err_perq_ci": analysis.wilson_ci(e1 + e2, 2 * a),
            "raw_err_perq": (raw_e1 + raw_e2) / (2 * (tot or 1))}


def eval_fusion(count: dict, basis: str, m: int, n: int) -> dict:
    """融合制备:先后选 v=b2=0(投影+制备验证),再后选 stabilizer 奇偶。
    码位 cb=(b0,b1,b3,b4);逻辑读出公式与 eval_422 同(cb 代入)。
    raw_* 指"仅 v 后选、无 stabilizer 后选"(v 前的态是 Bell⊗Bell,逻辑值未定义,
    故 raw 的基线就是 v 后选)。"""
    tot = v_acc = acc = e1 = e2 = ej = raw_e1 = raw_e2 = 0
    for key, val in count.items():
        b = analysis.bits_of(key, 5)
        tot += val
        if b[2] != 0:                       # 验证比特后选(融合投影 + 制备检错)
            continue
        v_acc += val
        cb = (b[0], b[1], b[3], b[4])
        if basis == "Z":
            l1, l2 = cb[0] ^ cb[2], cb[0] ^ cb[1]
        else:
            l1, l2 = cb[0] ^ cb[1], cb[0] ^ cb[2]
        raw_e1 += val * (l1 != m)
        raw_e2 += val * (l2 != n)
        if (cb[0] ^ cb[1] ^ cb[2] ^ cb[3]) == 0:  # stabilizer 后选
            acc += val
            e1 += val * (l1 != m)
            e2 += val * (l2 != n)
            ej += val * ((l1 != m) or (l2 != n))
    a = acc or 1
    return {"shots": tot,
            "v_acceptance": v_acc / (tot or 1),
            "acceptance": acc / (tot or 1),           # 总接受率(v ∧ stabilizer)
            "acceptance_ci": analysis.wilson_ci(acc, tot or 1),
            "ps_err_L1": e1 / a, "ps_err_L2": e2 / a, "ps_err_joint": ej / a,
            "ps_err_perq": (e1 + e2) / (2 * a),
            "ps_err_perq_ci": analysis.wilson_ci(e1 + e2, 2 * a),
            "raw_err_perq": (raw_e1 + raw_e2) / (2 * (v_acc or 1))}


def main():
    st = client.status()
    print(f"后端状态: {st}")
    ch = chipmod.load(CHIP, refresh=True)
    eta = load_today_eta()
    stars = pick_stars(ch, eta, N_STARS)
    paths = pick_paths(ch, eta, N_PATHS)
    # 基线扩到 5 个候选(当日 η top-5),报告时以"实测最好者"为最严苛基线,
    # 防止 η 排名与实测错误率不一致导致基线偏弱
    best_phys = sorted(eta, key=lambda q: -eta[q])[:5]
    print(f"选出 {len(stars)} 个星形:")
    for s in stars:
        print(f"  center={s['center']} leaves={s['leaves']} "
              f"min_eta={s['score_min']:.3f} cz={s['cz']}")
    print(f"选出 {len(paths)} 条融合 5-path:")
    for p in paths:
        print(f"  path={p['path']} v={p['vq']} min_eta={p['score_min']:.3f}")
    print(f"物理基线比特(η最高): {best_phys} eta={[round(eta[q],3) for q in best_phys]}")

    specs = []
    for s in stars:
        tq = [s["center"], *s["leaves"]]
        sid = f"c{s['center']}"
        for basis, builder in (("Z", circ_422_z), ("X", circ_422_x)):
            for m in (0, 1):
                for n in (0, 1):
                    specs.append(dict(
                        circuit=builder(m, n), chip=CHIP, shots=SHOTS,
                        name=f"c422_{sid}_{basis}{m}{n}", exp=EXP,
                        label=f"{sid}_{basis}_{m}{n}",
                        options={"target_qubits": tq},
                        meta={"kind": "code422", "star": s, "basis": basis,
                              "m": m, "n": n}))
        for basis in ("Z", "X"):          # 同substrate裸对照:星心比特裸存(与该星编码同批比特,相邻交错)
            for bit in (0, 1):
                specs.append(dict(
                    circuit=circ_phys(basis, bit), chip=CHIP, shots=SHOTS,
                    name=f"physss_{sid}_{basis}{bit}", exp=EXP,
                    label=f"physss_{sid}_{basis}_{bit}",
                    options={"target_qubits": [s["center"]]},
                    meta={"kind": "phys_ss", "qubit": s["center"], "star": s["center"],
                          "basis": basis, "bit": bit}))
    for p in paths:                       # 融合制备臂(FT 验证),与星形/物理同批交错
        pid = f"v{p['vq']}"
        for basis, builder in (("Z", circ_fusion_z), ("X", circ_fusion_x)):
            for m in (0, 1):
                for n in (0, 1):
                    specs.append(dict(
                        circuit=builder(m, n), chip=CHIP, shots=SHOTS,
                        name=f"c422f_{pid}_{basis}{m}{n}", exp=EXP,
                        label=f"{pid}_{basis}_{m}{n}",
                        options={"target_qubits": list(p["path"])},
                        meta={"kind": "fusion", "path": p, "basis": basis,
                              "m": m, "n": n}))
    for q in best_phys:
        for basis in ("Z", "X"):
            for bit in (0, 1):
                specs.append(dict(
                    circuit=circ_phys(basis, bit), chip=CHIP, shots=SHOTS,
                    name=f"phys_q{q}_{basis}{bit}", exp=EXP,
                    label=f"phys_q{q}_{basis}_{bit}",
                    options={"target_qubits": [q]},
                    meta={"kind": "phys", "qubit": q, "basis": basis, "bit": bit}))

    print(f"共 {len(specs)} 条线路 × {SHOTS} shots")
    recs = client.run_batch(specs, timeout=1200)

    # ---------------- 汇总
    rows = []
    for r in recs:
        c = client.counts(r)
        if not c:
            continue
        meta = r["meta"]
        if meta["kind"] == "code422":
            ev = eval_422(c, meta["basis"], meta["m"], meta["n"])
            ev = {kk: (round(v, 5) if isinstance(v, float)
                       else [round(x, 5) for x in v] if isinstance(v, (tuple, list))
                       else v) for kk, v in ev.items()}
            st = meta["star"]
            rows.append({**{k: meta[k] for k in ("basis", "m", "n")},
                         "star": st["center"],
                         "star_qubits": [st["center"], *st["leaves"]],
                         "label": r["label"], **ev})
        elif meta["kind"] == "fusion":
            ev = eval_fusion(c, meta["basis"], meta["m"], meta["n"])
            ev = {kk: (round(v, 5) if isinstance(v, float)
                       else [round(x, 5) for x in v] if isinstance(v, (tuple, list))
                       else v) for kk, v in ev.items()}
            p = meta["path"]
            rows.append({"kind": "fusion",
                         **{k: meta[k] for k in ("basis", "m", "n")},
                         "vq": p["vq"], "path_qubits": list(p["path"]),
                         "label": r["label"], **ev})
        else:  # phys(最强物理基线)或 phys_ss(同substrate,星心裸存)
            n1 = analysis.bit_p1(c, 1)[0]
            err = n1 if meta["bit"] == 0 else 1 - n1
            sh = analysis.shots_of(c)
            rows.append({"kind": meta["kind"], "qubit": meta["qubit"],
                         "star": meta.get("star"),
                         "basis": meta["basis"], "bit": meta["bit"],
                         "label": r["label"], "err": round(err, 5),
                         "err_ci": [round(x, 5)
                                    for x in analysis.wilson_ci(int(round(err * sh)), sh or 1)],
                         "shots": sh})

    OUT.mkdir(parents=True, exist_ok=True)
    summary = {"provenance": provenance.stamp(
                   anchor=None, chip=CHIP, tids=provenance.tids_of(recs),
                   qubits={"stars": [[s["center"], *s["leaves"]] for s in stars],
                           "fusion_paths": [list(p["path"]) for p in paths],
                           "phys_baseline": best_phys},
                   extra={"anchor": "e14_selfcal_eta + pick_stars/pick_paths",
                          "chip_calibration_time": ch.calibration_time,
                          "measured_date": TODAY}),
               "rows": rows}
    (OUT / f"summary_{RUN}.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2))

    print("\n===== [[4,2,2]] 对决(后选逻辑错误率/每逻辑比特)=====")
    stars_seen = sorted({r["star"] for r in rows
                         if "ps_err_perq" in r and "star" in r})  # 只 code422 行
    for sid in stars_seen:
        sub = [r for r in rows if r.get("star") == sid and "ps_err_perq" in r]
        for basis in ("Z", "X"):
            bs = [r for r in sub if r["basis"] == basis]
            if not bs:
                continue
            acc = sum(r["acceptance"] for r in bs) / len(bs)
            ps = sum(r["ps_err_perq"] for r in bs) / len(bs)
            raw = sum(r["raw_err_perq"] for r in bs) / len(bs)
            print(f"  star c{sid} {basis}基: raw={raw:.4f} → 后选={ps:.4f} "
                  f"(acceptance={acc:.3f})")
    fusion_seen = sorted({r["vq"] for r in rows if r.get("kind") == "fusion"})
    for vq in fusion_seen:
        sub = [r for r in rows if r.get("kind") == "fusion" and r["vq"] == vq]
        for basis in ("Z", "X"):
            bs = [r for r in sub if r["basis"] == basis]
            if not bs:
                continue
            acc = sum(r["acceptance"] for r in bs) / len(bs)
            vac = sum(r["v_acceptance"] for r in bs) / len(bs)
            ps = sum(r["ps_err_perq"] for r in bs) / len(bs)
            raw = sum(r["raw_err_perq"] for r in bs) / len(bs)
            print(f"  fusion v{vq} {basis}基: raw(v后)={raw:.4f} → 后选={ps:.4f} "
                  f"(v_acc={vac:.3f}, 总acc={acc:.3f})")
    phys = [r for r in rows if r.get("kind") == "phys"]
    for basis in ("Z", "X"):
        bs = [r for r in phys if r["basis"] == basis]
        if bs:
            e = sum(r["err"] for r in bs) / len(bs)
            print(f"  最强物理基线 {basis}基: err={e:.4f} ({[r['qubit'] for r in bs]})")
    physss = [r for r in rows if r.get("kind") == "phys_ss"]
    for basis in ("Z", "X"):
        bs = [r for r in physss if r["basis"] == basis]
        if bs:
            e = sum(r["err"] for r in bs) / len(bs)
            print(f"  同substrate裸(星心) {basis}基: err={e:.4f}")
    print(f"\n已写 {OUT / f'summary_{RUN}.json'}")


if __name__ == "__main__":
    main()
