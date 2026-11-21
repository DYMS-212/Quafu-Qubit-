"""e25 · Iceberg 检错码族 [[k+2, k, 2]] 规模扫描:编码率 → 1 时买到什么、付出什么。

动机:[[4,2,2]](e10)是本码族 k=2 的特例。同一族码把块做大(k=2→8,
n=4→10),编码率 k/n 从 1/2 升到 4/5——**"错误率 × 接受率 × 编码率"的三方
权衡曲线**是真机上罕见的规模化叙事,且平台约束(末端测量+后选)完美兼容。

构造(全部 statevector 自校准,绝不手推):
  稳定子 S_X=X⊗ⁿ, S_Z=Z⊗ⁿ;逻辑 X̄ᵢ=X₀Xᵢ, Z̄ᵢ=ZᵢZ_{n-1}(i=1..k)。
  |0̄…0̄⟩ = GHZ_n(路径中点向两侧 CX,深度 ⌈(n-1)/2⌉);|1̄…1̄⟩ = 再打一层
  物理 X(Πᵢ X̄ᵢ 的支撑)。存储期 K 层远端 CZ 负载(同 e11 机制,真 idle)。
  末端全 Z 测量:S_Z 奇偶后选(单比特翻转可检)→ 逻辑读出 z̄ᵢ = bᵢ⊕b_{n-1}。
  口径:bit-flip 检错(与 repetition/e10-Z 同 caliber);相位保护的完整口径
  由 [[4,2,2]](e10/e24)承担,本实验专攻"率-开销"轴。

对照:同链前 k 个物理比特裸存(同 K 负载、同批交错)。
产出:每 (k, L, K) 一行:acceptance / 每逻辑比特错误率 / 任一逻辑错 / 编码率。

运行: ./qrun experiments/e25_iceberg.py [chip]      # 默认 Baihua
      ./qrun experiments/e25_iceberg.py selftest    # 零额度 statevector 自检
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, anchors as anchmod, client, provenance  # noqa: E402
from qfab.qasm import QASM  # noqa: E402

CHIP = sys.argv[1] if len(sys.argv) > 1 else "Baihua"
TODAY = time.strftime("%Y-%m-%d")
RUN = time.strftime("%Y-%m-%d_%H%M")
TAG = RUN if CHIP == "Baihua" else f"{CHIP}_{RUN}"
EXP = f"e25_iceberg/{TAG}"
OUT = Path(__file__).resolve().parent.parent / "results" / "e25_iceberg"
SHOTS = 4096
KS = [2, 4, 6, 8]           # 逻辑比特数;n = k+2 ∈ {4,6,8,10}
K_IDLE = [0, 8]             # 负载层数(0=静态存读,8=真 idle 窗口)


# ------------------------------------------------------------------ 线路
def xbar_pattern(k: int) -> list[int]:
    """Πᵢ₌₁..k X̄ᵢ 的物理支撑:X₀^k · X₁…X_k(mod 2)。"""
    supp = list(range(1, k + 1))
    if k % 2 == 1:
        supp = [0] + supp
    return supp


def circ(k: int, logical: int, K: int) -> tuple[str, int]:
    """[[k+2,k,2]] 线路。返回 (qasm, 总宽度)。logical∈{0,1}=全0̄/全1̄。
    GHZ 从路径中点向两侧生长(深度减半);K>0 时挂 2 个负载比特。"""
    n = k + 2
    width = n + (2 if K > 0 else 0)
    qc = QASM(width)
    mid = n // 2
    qc.h(mid)
    for i in range(mid, n - 1):          # 右臂
        qc.cx(i, i + 1)
    for i in range(mid, 0, -1):          # 左臂
        qc.cx(i, i - 1)
    if logical:
        for q in xbar_pattern(k):
            qc.x(q)
    qc.barrier()
    if K > 0:
        la, lb = n, n + 1
        for _ in range(K):
            qc.cz(la, lb)
            qc.barrier()
    qc.measure_all()
    return qc.qasm(), width


def circ_bare(k: int, logical: int, K: int) -> tuple[str, int]:
    """裸对照:k 个物理比特同存同测(同 K 负载)。"""
    width = k + (2 if K > 0 else 0)
    qc = QASM(width)
    if logical:
        for q in range(k):
            qc.x(q)
    qc.barrier()
    if K > 0:
        for _ in range(K):
            qc.cz(k, k + 1)
            qc.barrier()
    qc.measure_all()
    return qc.qasm(), width


# ------------------------------------------------------------------ 分析
def extract(count: dict, k: int, K: int, exp_parity: int,
            exp_logical: list[int]) -> dict:
    """S_Z 奇偶后选 → 逻辑读出。返回 acceptance / 每逻辑错误率 / 任一逻辑错。"""
    n = k + 2
    width = n + (2 if K > 0 else 0)
    tot = acc = any_err = 0
    per = [0] * k
    for key, v in count.items():
        b = analysis.bits_of(key, width)
        tot += v
        if (sum(b[:n]) % 2) != exp_parity:
            continue
        acc += v
        bad = False
        for i in range(1, k + 1):
            if (b[i] ^ b[n - 1]) != exp_logical[i - 1]:
                per[i - 1] += v
                bad = True
        if bad:
            any_err += v
    per_mean = sum(per) / (k * acc) if acc else None
    return {"acceptance": round(acc / (tot or 1), 4),
            "n_acc": acc,
            "err_per_logical": round(per_mean, 6) if acc else None,
            "err_any": round(any_err / acc, 6) if acc else None,
            "per_logical": [round(x / acc, 6) if acc else None for x in per]}


def ideal(k: int, logical: int) -> tuple[int, list[int]]:
    """statevector 自校准:理想 S_Z 奇偶 + 理想逻辑读出(与线路永远同步)。"""
    from tools.verify_circuits_offline import simulate
    cnt = simulate(circ(k, logical, 0)[0])
    n = k + 2
    parities = set()
    logicals = set()
    for key in cnt:
        b = analysis.bits_of(key, n)
        parities.add(sum(b[:n]) % 2)
        logicals.add(tuple(b[i] ^ b[n - 1] for i in range(1, k + 1)))
    assert len(parities) == 1, f"k={k} L{logical}: S_Z 奇偶不定 {parities}"
    assert len(logicals) == 1, f"k={k} L{logical}: 逻辑读出不定 {logicals}"
    return parities.pop(), list(logicals.pop())


def selftest():
    for k in KS:
        for logical in (0, 1):
            p, lg = ideal(k, logical)
            exp = [logical] * k
            ok = lg == exp
            print(f"  {'✅' if ok else '❌'} [[{k + 2},{k},2]] L{logical}: "
                  f"S_Z奇偶={p} 逻辑读出={lg}")
            if not ok:
                sys.exit(1)
            from tools.verify_circuits_offline import simulate
            m = extract(simulate(circ(k, logical, 0)[0]), k, 0, p, lg)
            assert m["acceptance"] == 1.0 and m["err_any"] == 0.0, m
    print("selftest 全部通过 ✅")


# ------------------------------------------------------------------ 主流程
def main():
    if CHIP == "selftest":
        selftest()
        return
    st = client.status()
    print(f"后端状态: {st}")
    if not isinstance(st.get(CHIP), int):
        print(f"⚠️ {CHIP} 当前不可用({st.get(CHIP)}),中止。")
        return
    A = anchmod.latest(CHIP)
    chains = A.get("chains") or {}
    chain = max(chains.values(), key=len, default=[])
    load = list(A.get("load_edge") or [])
    n_max = max(k + 2 for k in KS)
    if len(chain) < n_max or len(load) != 2:
        print(f"⚠️ 锚点链不足 {n_max} 比特或缺负载边,可跑的 k 将截断。")

    cal = {}          # (k, L) → (parity, logical)
    specs = []
    for k in KS:
        n = k + 2
        if len(chain) < n:
            continue
        for logical in (0, 1):
            cal[(k, logical)] = ideal(k, logical)
            for K in K_IDLE:
                q, w = circ(k, logical, K)
                tq = chain[:n] + (load if K > 0 else [])
                specs.append(dict(circuit=q, chip=CHIP, shots=SHOTS,
                                  name=f"ice_k{k}L{logical}K{K}", exp=EXP,
                                  label=f"ice_k{k}_L{logical}_K{K}",
                                  options={"target_qubits": tq},
                                  meta={"kind": "ice", "k": k, "logical": logical,
                                        "K": K}))
                qb, wb = circ_bare(k, logical, K)
                tqb = chain[:k] + (load if K > 0 else [])
                specs.append(dict(circuit=qb, chip=CHIP, shots=SHOTS,
                                  name=f"iceb_k{k}L{logical}K{K}", exp=EXP,
                                  label=f"bare_k{k}_L{logical}_K{K}",
                                  options={"target_qubits": tqb},
                                  meta={"kind": "bare", "k": k, "logical": logical,
                                        "K": K}))
    print(f"共 {len(specs)} 条线路 × {SHOTS} shots(链 {chain},负载 {load})")
    recs = client.run_batch(specs, timeout=1200)

    rows = []
    for r in recs:
        c = client.counts(r)
        if not c:
            continue
        m = r["meta"]
        k = m["k"]
        if m["kind"] == "ice":
            p, lg = cal[(k, m["logical"])]
            out = extract(c, k, m["K"], p, lg)
            if out["n_acc"]:
                lo, hi = analysis.wilson_ci(
                    round((out["err_per_logical"] or 0) * k * out["n_acc"]),
                    k * out["n_acc"])
                out["err_ci"] = [round(lo, 6), round(hi, 6)]
            rows.append({**m, "chip": CHIP, "n": k + 2,
                         "rate": round(k / (k + 2), 3),
                         "qubits": chain[:k + 2], "tid": r.get("tid"), **out})
        else:
            width = k + (2 if m["K"] > 0 else 0)
            tot = analysis.shots_of(c)
            errs = 0
            for key, v in c.items():
                b = analysis.bits_of(key, width)
                errs += v * sum(1 for i in range(k) if b[i] != m["logical"])
            per = errs / (k * tot)
            lo, hi = analysis.wilson_ci(round(errs), k * tot)
            rows.append({**m, "chip": CHIP, "qubits": chain[:k],
                         "tid": r.get("tid"),
                         "err_per_logical": round(per, 6),
                         "err_ci": [round(lo, 6), round(hi, 6)]})

    OUT.mkdir(parents=True, exist_ok=True)
    summary = {"provenance": provenance.stamp(
                   anchor=A, chip=CHIP, tids=provenance.tids_of(recs),
                   qubits={"chain": chain, "load_edge": load},
                   extra={"caliber": "bit-flip 检错口径(S_Z 奇偶后选);"
                                     "相位口径由 [[4,2,2]] e10/e24 承担",
                          "err_ci_note": "per-logical 汇总按 k×n_acc 伯努利近似"
                                         "(逻辑位间相关,CI 偏窄,判定以跨日为准)"}),
               "rows": rows}
    (OUT / f"summary_{TAG}.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2))

    print("\n===== Iceberg [[k+2,k,2]] 规模扫描 =====")
    print(f"{'k':>3s} {'率':>6s} {'K':>3s} {'L':>2s} {'接受率':>8s} "
          f"{'每逻辑错':>10s} {'裸对照':>10s}")
    for k in KS:
        for K in K_IDLE:
            for logical in (0, 1):
                ice = [r for r in rows if r["kind"] == "ice" and r["k"] == k
                       and r["K"] == K and r["logical"] == logical]
                bare = [r for r in rows if r["kind"] == "bare" and r["k"] == k
                        and r["K"] == K and r["logical"] == logical]
                if ice:
                    i0 = ice[0]
                    b0 = bare[0] if bare else {}
                    print(f"{k:>3d} {i0['rate']:>6.2f} {K:>3d} {logical:>2d} "
                          f"{i0['acceptance']:>8.3f} "
                          f"{i0['err_per_logical'] if i0['err_per_logical'] is not None else float('nan'):>10.5f} "
                          f"{b0.get('err_per_logical', float('nan')):>10.5f}")
    print(f"已写 {OUT / f'summary_{TAG}.json'}")


if __name__ == "__main__":
    main()
