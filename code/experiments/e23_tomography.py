"""e23 · 纠缠态断层成像:Bell(2 比特,9 基)+ GHZ(3 比特,27 基)完整密度矩阵。

动机:密度矩阵"城市图"(Re ρ 四柱 + Im ρ 平地)是最直观的"纠缠真实存在"证据,
观众一眼懂;同一套层析跑三台机器 = 跨机对比的门面图;raw vs 读出缓解双版本
= 缓解板块的可视化。

方法:制备(Bell: H+CX;GHZ3: H+CX+CX 链)→ 每比特测量基旋转(X:H;Y:S†H;Z:无)
→ 全部 3^n 个基组合 × 2048 shots → Pauli 期望线性反演 ρ → 最近物理态投影
(PSD 截断+归一)。指标:保真度 F=⟨ψ|ρ|ψ⟩、纯度 Tr ρ²、负性(2q,纠缠的数学判据)。
有当日自校准的机器同时给读出缓解版(张量反演,qfab.mitigation)。

选址:Baihua=当日 atlas 最优边 / 当日锚点星;Dongling/Shenglian=其 atlas 或
厂商最优(三机全部定点放置;Shenglian 标定经 chip-details 端点,见 qfab/client.py)。

运行: ./qrun experiments/e23_tomography.py [chip]      # 默认 Baihua
      ./qrun experiments/e23_tomography.py selftest    # 零额度,statevector 自检
"""
from __future__ import annotations

import json
import sys
import time
from itertools import product
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, client, mitigation, provenance  # noqa: E402
from qfab.qasm import QASM  # noqa: E402

PROJ = Path(__file__).resolve().parent.parent
CHIP = sys.argv[1] if len(sys.argv) > 1 else "Baihua"
RUN = time.strftime("%Y-%m-%d_%H%M")
TODAY = time.strftime("%Y-%m-%d")
EXP = f"e23_tomography/{RUN}"
OUT = PROJ / "results" / "e23_tomography"
SHOTS = 2048

_P = {"I": np.eye(2), "X": np.array([[0, 1], [1, 0]], dtype=complex),
      "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.diag([1, -1]).astype(complex)}


# ------------------------------------------------------------------ 线路
def tomo_circ(n: int, setting: str) -> str:
    """制备 n 比特 GHZ 型纠缠态(n=2 即 Bell)+ 按 setting 旋转测量基。
    setting[i] = 比特 i 的测量基(X/Y/Z)。"""
    qc = QASM(n)
    qc.h(0)
    for i in range(n - 1):
        qc.cx(i, i + 1)
    qc.barrier()
    for i, b in enumerate(setting):
        if b == "X":
            qc.h(i)
        elif b == "Y":
            qc.sdg(i).h(i)
    qc.barrier().measure_all()
    return qc.qasm()


# ------------------------------------------------------------------ 重建
def reconstruct(counts_by_setting: dict[str, dict], n: int) -> np.ndarray:
    """Pauli 期望线性反演。张量序:kron(比特 n-1, ..., 比特 0)
    ↔ 基态索引 = b_{n-1}...b_0(与平台 count 位序一致)。"""
    rho = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for ps in product("IXYZ", repeat=n):          # ps[i] = 比特 i 的 Pauli
        if all(p == "I" for p in ps):
            e = 1.0
        else:
            vals = []
            for setting, cnt in counts_by_setting.items():
                if not cnt:
                    continue
                if all(p == "I" or p == s for p, s in zip(ps, setting)):
                    tot = sum(cnt.values())
                    acc = 0.0
                    for key, v in cnt.items():
                        bits = analysis.bits_of(key, n)
                        sign = (-1) ** sum(bits[i] for i, p in enumerate(ps)
                                           if p != "I")
                        acc += sign * v
                    vals.append(acc / tot)
            if not vals:
                continue
            e = float(np.mean(vals))
        op = _P[ps[n - 1]]
        for i in range(n - 2, -1, -1):
            op = np.kron(op, _P[ps[i]])
        rho += e * op
    return rho / (2 ** n)


def project_psd(rho: np.ndarray) -> np.ndarray:
    """最近物理态投影:特征值截负 + 归一。"""
    w, v = np.linalg.eigh((rho + rho.conj().T) / 2)
    w = np.clip(w, 0, None)
    w = w / w.sum()
    return (v * w) @ v.conj().T


def _pauli_exp(rho: np.ndarray, ps: str) -> float:
    n = len(ps)
    op = _P[ps[n - 1]]
    for i in range(n - 2, -1, -1):
        op = np.kron(op, _P[ps[i]])
    return float(np.real(np.trace(rho @ op)))


def metrics(rho: np.ndarray, n: int) -> dict:
    psi = np.zeros(2 ** n)
    psi[0] = psi[-1] = 1 / np.sqrt(2)             # (|0..0>+|1..1>)/√2
    out = {"fidelity": float(np.real(psi @ rho @ psi)),
           "purity": float(np.real(np.trace(rho @ rho)))}
    if n == 2:                                    # 负性(部分转置最负特征值)
        pt = rho.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
        out["negativity"] = float(max(0.0, -np.linalg.eigvalsh(pt).min()))
        # CHSH 最大违背(Horodecki 判据):S_max=2√(t1+t2),经典界 2,量子上限 2√2
        T = np.array([[_pauli_exp(rho, a + b) for b in "XYZ"] for a in "XYZ"])
        ev = sorted(np.linalg.eigvalsh(T.T @ T))[::-1]
        out["chsh_max"] = float(2 * np.sqrt(max(0.0, ev[0] + ev[1])))
    if n == 3:                                    # Mermin 算符,经典界 2、量子上限 4
        out["mermin"] = float(_pauli_exp(rho, "XXX") - _pauli_exp(rho, "XYY")
                              - _pauli_exp(rho, "YXY") - _pauli_exp(rho, "YYX"))
    return out


# ------------------------------------------------------------------ 选址
def pick_sites():
    """返回 (pair, triple):自测图谱(e03)优先,否则厂商标定最优。"""
    from qfab import anchors as anchmod, chip as chipmod
    atlas = sorted((PROJ / "results" / "e03_bell").glob(f"atlas_{CHIP}_*.json"))
    pair = None
    if atlas:
        rows = json.loads(atlas[-1].read_text()).get("rows", [])
        best = max((r for r in rows if r.get("bell_mean")),
                   key=lambda r: r["bell_mean"], default=None)
        if best:
            pair = list(best["edge"])
    if CHIP == "Baihua":
        a = anchmod.latest(CHIP)
        st = a.get("star") or {}
        triple = [st["leaves"][0], st["center"], st["leaves"][1]] \
            if st.get("leaves") else None
        return pair or [67, 68], triple
    ch = chipmod.load(CHIP, refresh=True)         # 其他有拓扑机:厂商最优
    av = set(ch.available_qubits())
    edges = [(a, b) for a in av for b in ch.neighbors(a)
             if b in av and a < b and (ch.edge_fid(a, b) or 0) > 0.5]
    if not pair and edges:
        pair = list(max(edges, key=lambda e: ch.edge_fid(*e)))
    triple = None
    for a, b in edges:                             # 链 a-b-c
        for c in ch.neighbors(b):
            if c in av and c != a and (ch.edge_fid(b, c) or 0) > 0.5:
                cand = [a, b, c]
                s = min(ch.edge_fid(a, b), ch.edge_fid(b, c))
                if triple is None or s > triple[0]:
                    triple = (s, cand)
    return pair, (triple[1] if triple else None)


def cal_for(phys: list[int]) -> list[tuple[float, float]] | None:
    """当日自校准 (F0,F1) 列表(读出缓解用);无数据 → None。"""
    try:
        if CHIP == "Baihua":
            f = sorted((PROJ / "results" / "e02_selfcal")
                       .glob(f"assignment_{TODAY}*.json"))[-1]
            cal = {r["qubit"]: r for r in json.loads(f.read_text())}
        else:
            d = sorted((PROJ / "results" / "e15_crossbackend")
                       .glob(f"{CHIP}_*"))[-1]
            pairs: dict[tuple, dict[int, dict]] = {}
            for f in sorted(d.glob("selfcal_prep*__*.json")):
                r = json.loads(f.read_text())
                if not r.get("result", {}).get("count"):
                    continue
                p = 1 if r["meta"]["prep"] == "prep1" else 0
                pairs.setdefault(tuple(r["meta"]["qubits"]), {})[p] = r
            cal = {}
            for qs, pr in pairs.items():
                if 0 not in pr or 1 not in pr:
                    continue
                asg = analysis.assignment_from_p1(
                    analysis.bit_p1(pr[0]["result"]["count"], len(qs)),
                    analysis.bit_p1(pr[1]["result"]["count"], len(qs)))
                cal.update({q: asg[i] for i, q in enumerate(qs)})
        return [(cal[q]["F0"], cal[q]["F1"]) for q in phys]
    except Exception:  # noqa: BLE001
        return None


# ------------------------------------------------------------------ 主流程
def analyze(name: str, n: int, recs_by_setting: dict[str, dict],
            phys: list[int] | None) -> dict:
    rho = project_psd(reconstruct(recs_by_setting, n))
    out = {"n": n, "qubits": phys, "settings": len(recs_by_setting),
           "shots": SHOTS,
           "rho_re": np.round(rho.real, 5).tolist(),
           "rho_im": np.round(rho.imag, 5).tolist(), **metrics(rho, n)}
    if phys:
        fs = cal_for(phys)
        if fs:
            mit = {s: mitigation.correct_counts(c, fs)[0]
                   for s, c in recs_by_setting.items() if c}
            rho_m = project_psd(reconstruct(mit, n))
            out["mitigated"] = {"rho_re": np.round(rho_m.real, 5).tolist(),
                                "rho_im": np.round(rho_m.imag, 5).tolist(),
                                **metrics(rho_m, n)}
    return out


def selftest():
    from tools.verify_circuits_offline import simulate
    for n in (2, 3):
        cbs = {"".join(s): simulate(tomo_circ(n, "".join(s)))
               for s in product("XYZ", repeat=n)}
        rho = project_psd(reconstruct(cbs, n))
        m = metrics(rho, n)
        ok = m["fidelity"] > 0.999 and m["purity"] > 0.999
        if n == 2:
            ok &= abs(m["chsh_max"] - 2 * np.sqrt(2)) < 1e-6
            extra = f"negativity={m['negativity']:.4f} CHSH={m['chsh_max']:.4f}"
        else:
            ok &= abs(abs(m["mermin"]) - 4) < 1e-6
            extra = f"Mermin={m['mermin']:.4f}"
        print(f"  {'✅' if ok else '❌'} n={n}: F={m['fidelity']:.6f} "
              f"purity={m['purity']:.6f} {extra}")
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
    pair, triple = pick_sites()
    print(f"选址: bell={pair} ghz3={triple}(None=跳过)")

    specs = []
    plan = [("bell", 2, pair)]
    if triple:
        plan.append(("ghz3", 3, triple))
    for name, n, phys in plan:
        for s in product("XYZ", repeat=n):
            setting = "".join(s)
            spec = dict(circuit=tomo_circ(n, setting), chip=CHIP, shots=SHOTS,
                        name=f"tomo_{name}_{setting}", exp=EXP,
                        label=f"{name}_{setting}",
                        meta={"state": name, "n": n, "setting": setting,
                              "phys": phys})
            if phys:
                spec["options"] = {"target_qubits": phys}
            specs.append(spec)
    print(f"共 {len(specs)} 条线路 × {SHOTS} shots")
    recs = client.run_batch(specs, timeout=1200)

    states = {}
    for name, n, phys in plan:
        cbs = {}
        for r in recs:
            if r["meta"]["state"] != name or not client.counts(r):
                continue
            cbs[r["meta"]["setting"]] = client.counts(r)
        if len(cbs) < 3 ** n:
            print(f"⚠️ {name}: 仅 {len(cbs)}/{3**n} 个基有数据")
        if not cbs:
            continue
        states[name] = analyze(name, n, cbs, phys)
        m = states[name]
        print(f"  {name}: F={m['fidelity']:.4f} purity={m['purity']:.4f}"
              + (f" negativity={m['negativity']:.4f}" if "negativity" in m else "")
              + (f" | 缓解后 F={m['mitigated']['fidelity']:.4f}"
                 if "mitigated" in m else ""))

    OUT.mkdir(parents=True, exist_ok=True)
    summary = {"provenance": provenance.stamp(
                   anchor=None, chip=CHIP, tids=provenance.tids_of(recs),
                   qubits={"bell": pair, "ghz3": triple},
                   extra={"mode": "targeted"}),
               "chip": CHIP, "states": states}
    (OUT / f"summary_{CHIP}_{RUN}.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"已写 {OUT / f'summary_{CHIP}_{RUN}.json'}")


if __name__ == "__main__":
    main()
