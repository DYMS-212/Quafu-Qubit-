"""e11 · 编码存储寿命:p_L(d, t_idle) 与逻辑寿命增益(headline)。

动机:粉碎"repetition=读出平均"的质疑——idle 存储期发生的是**真实退相干错误**
(T1 衰减→bit flip),若编码在扫描存储时长下仍压制逻辑错误,即证明是在纠
"存储中发生的错误",不是平均读出噪声。

方法(无 mid-circuit 平台上最接近 memory 实验的诚实版本):
  编码(bit 码 d=3/5/7,链锚点由当日选址给出)→ K 层 CZ 负载造真 idle
  (负载打在远处边,每层 barrier 全域同步)
  → syndrome 提取 → 末端测量 → 多数表决(解码器阶梯 e12 已证 T=1 最优≈表决)。
  物理基线:同 K 负载下单比特裸存(链上 q8 + 全芯片最好 q150,|0>/|1>)。

**phase 臂**:bit 码存 |0>/|1> 保护的是经典信息(T1 翻转);
为堵"经典冗余"质疑,增加 phase-flip 码存 |+>/|−>(数据比特 H 夹层,存储期
积累的是 T2 退相位),同 K 负载、同批交错、同一套多数表决(H 后回 Z 基)。
物理基线同样加 |+>/|−>(Ramsey 式裸存)。若相关退相位破坏错误独立性导致
编码失效,如实报告——那是"相位存储的边界"结果。

产出指标:
  - p_L(d, K) 曲线族 + 物理 p(K) → 逻辑比特寿命增益(达到 10% 错误的时间比)
  - Λ(d→d+2, K) = p_L(d)/p_L(d+2) 随存储时长的衰减 → 距离增益在哪失效
  - L0/L1 分开报告(T1 偏置噪声,诚实处理);bit/phase 码分开报告

运行: ./qrun experiments/e11_encoded_memory.py [chip] [run|analyze|fixk0]
     # 默认 Baihua/run;非 Baihua 机 = lite 档(K∈{0,8,32},其余同口径)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, anchors as anchmod, client, repcode, provenance  # noqa: E402
from qfab.qasm import QASM  # noqa: E402

_CHIPS = ("Baihua", "Dongling", "Shenglian")
CHIP = next((a for a in sys.argv[1:] if a in _CHIPS), "Baihua")
MODE = next((a for a in sys.argv[1:] if a in ("run", "analyze", "fixk0")), "run")
TODAY = time.strftime("%Y-%m-%d")
RUN = time.strftime("%Y-%m-%d_%H%M")   # 同日多轮互不覆盖
TAG = RUN if CHIP == "Baihua" else f"{CHIP}_{RUN}"
EXP = f"e11_memory/{TAG}"
OUT = Path(__file__).resolve().parent.parent / "results" / "e11_memory"
SHOTS = 4096

_A = anchmod.latest(CHIP)
CHAINS = {int(k): v for k, v in _A["chains"].items() if int(k) <= 7}  # 存储对决 d≤7 足够
LOAD_EDGE = tuple(_A["load_edge"])
# CZ 负载层数(≈存储时长);非 Baihua 机 lite 档三点即可拉出曲线
K_GRID = [0, 4, 8, 16, 32] if CHIP == "Baihua" else [0, 8, 32]
_cand = [CHAINS[3][0]] + list(_A["phys_baseline"])
PHYS_QUBITS = list(dict.fromkeys(q for q in _cand if q not in LOAD_EDGE))[:2]  # 链代表+最好,且不撞负载边


def build_encoded(d: int, bit: int, K: int, phase: bool = False) -> tuple[str, list[int]]:
    """repetition 链 + 2 负载比特(逻辑位 n-2,n-1)。phase=True:数据比特 H 夹层,
    存储期在 X 基(积累 T2 退相位),提取前转回 Z 基,解码仍是多数表决。"""
    n = 2 * d - 1
    qc = QASM(n + 2)
    la, lb = n, n + 1
    data = [2 * k for k in range(d)]
    if bit:
        for q in data:
            qc.x(q)
    if phase:
        for q in data:
            qc.h(q)
    qc.barrier()
    for _ in range(K):                     # 存储期:负载边打 CZ,其余比特真 idle
        qc.cz(la, lb)
        qc.barrier()
    if phase:                              # 转回 Z 基(barrier 已隔开,编译器不可对消)
        for q in data:
            qc.h(q)
    for k in range(d - 1):                 # syndrome 提取
        qc.cx(2 * k, 2 * k + 1)
        qc.cx(2 * (k + 1), 2 * k + 1)
    qc.barrier()
    qc.measure_all()
    return qc.qasm(), CHAINS[d] + list(LOAD_EDGE)


def build_phys(bit: int, K: int, q: int, phase: bool = False) -> tuple[str, list[int]]:
    if K == 0:
        # 坑:平台拒收"含完全无门且不测量的比特"的线路 → K=0 用纯 1 比特
        qc = QASM(1)
        if bit:
            qc.x(0)
        if phase:
            qc.h(0)
        qc.barrier()
        if phase:
            qc.h(0)
        qc.measure(0, 0)
        return qc.qasm(), [q]
    qc = QASM(3)                           # q0=被测,q1/q2=负载
    if bit:
        qc.x(0)
    if phase:
        qc.h(0)
    qc.barrier()
    for _ in range(K):
        qc.cz(1, 2)
        qc.barrier()
    if phase:
        qc.h(0)
    qc.measure(0, 0)
    return qc.qasm(), [q, *LOAD_EDGE]


def phys_specs(only_k0: bool = False):
    specs = []
    for q in PHYS_QUBITS:
        for phase in (False, True):
            code = "phase" if phase else "bit"
            tag = "P" if phase else ""
            for bit in (0, 1):
                for K in ([0] if only_k0 else K_GRID):
                    circ, tq = build_phys(bit, K, q, phase=phase)
                    specs.append(dict(circuit=circ, chip=CHIP, shots=SHOTS,
                                      name=f"memphys_q{q}{tag}L{bit}K{K}", exp=EXP,
                                      label=f"phys{q}{tag}_L{bit}_K{K}",
                                      options={"target_qubits": tq},
                                      meta={"kind": "phys", "qubit": q, "bit": bit,
                                            "K": K, "code": code}))
    return specs


def load_day_records() -> list[dict]:
    root = Path(__file__).resolve().parent.parent / "results" / "e11_memory"
    pref = "20" if CHIP == "Baihua" else CHIP
    days = sorted(d for d in root.iterdir()
                  if d.is_dir() and d.name.startswith(pref))
    recs = [json.loads(f.read_text()) for f in sorted(days[-1].glob("*.json"))]
    return [r for r in recs if (r.get("result") or {}).get("count")]


def main():
    mode = MODE
    if mode == "analyze":
        recs = load_day_records()
        print(f"离线分析 {len(recs)} 条已有记录")
    elif mode == "fixk0":
        print(f"后端状态: {client.status()}")
        client.run_batch(phys_specs(only_k0=True), timeout=1200)
        recs = load_day_records()
    else:
        print(f"后端状态: {client.status()}")
        specs = []
        for d in CHAINS:
            for phase in (False, True):
                code = "phase" if phase else "bit"
                tag = "P" if phase else ""
                for bit in (0, 1):
                    for K in K_GRID:
                        circ, tq = build_encoded(d, bit, K, phase=phase)
                        specs.append(dict(circuit=circ, chip=CHIP, shots=SHOTS,
                                          name=f"mem_d{d}{tag}L{bit}K{K}", exp=EXP,
                                          label=f"d{d}{tag}_L{bit}_K{K}",
                                          options={"target_qubits": tq},
                                          meta={"kind": "enc", "d": d, "bit": bit,
                                                "K": K, "code": code}))
        specs += phys_specs()
        print(f"共 {len(specs)} 条线路 × {SHOTS} shots")
        recs = client.run_batch(specs, timeout=1200)

    rows = []
    for r in recs:
        c = client.counts(r)
        if not c:
            continue
        m = r["meta"]
        if m["kind"] == "enc":
            d = m["d"]
            # 忽略最高 2 位(负载比特):把 count 串右侧 2d-1 位取出
            n = 2 * d - 1
            sub: dict = {}
            for key, v in c.items():
                kk = key.zfill(n + 2)[2:]
                sub[kk] = sub.get(kk, 0) + v
            err = repcode.logical_error_majority(sub, d, m["bit"])
            ci = analysis.wilson_ci(round(err * analysis.shots_of(sub)),
                                    analysis.shots_of(sub))
            rows.append({**m, "chip": CHIP, "qubits": CHAINS[d] + list(LOAD_EDGE),
                         "tid": r.get("tid"), "err": round(err, 5),
                         "ci": [round(x, 5) for x in ci]})
        else:
            p1 = analysis.bit_p1(c, 3)[0]
            err = p1 if m["bit"] == 0 else 1 - p1
            ci = analysis.wilson_ci(round(err * analysis.shots_of(c)),
                                    analysis.shots_of(c))
            rows.append({**m, "chip": CHIP, "qubits": ([m["qubit"], *LOAD_EDGE] if m["K"] else [m["qubit"]]),
                         "tid": r.get("tid"), "err": round(err, 5),
                         "ci": [round(x, 5) for x in ci]})

    OUT.mkdir(parents=True, exist_ok=True)
    summary = {"provenance": provenance.stamp(
                   anchor=_A, chip=CHIP, tids=provenance.tids_of(recs),
                   qubits={"chains": CHAINS, "load_edge": list(LOAD_EDGE),
                           "phys": PHYS_QUBITS}),
               "rows": rows}
    (OUT / f"summary_{TAG}.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2))

    for code in ("bit", "phase"):
        sub_rows = [r for r in rows if r.get("code", "bit") == code]
        if not sub_rows:
            continue
        desc = ("L1 = 存 |1>,T1 衰减主导" if code == "bit"
                else "L1 = 存 |−>,T2 退相位主导")
        print(f"\n===== p_L(d, K) 编码存储 · {code} 码({desc})=====")
        print(f"{'K':>4s}", *[f"{'d='+str(d):>10s}" for d in CHAINS],
              *[f"{'phys_q'+str(q):>10s}" for q in PHYS_QUBITS])
        for bit in (1, 0):
            print(f"-- L{bit} --")
            for K in K_GRID:
                vals = []
                for d in CHAINS:
                    v = [r for r in sub_rows if r["kind"] == "enc" and r["d"] == d
                         and r["bit"] == bit and r["K"] == K]
                    vals.append(f"{v[0]['err']:>10.4f}" if v else f"{'—':>10s}")
                for q in PHYS_QUBITS:
                    v = [r for r in sub_rows if r["kind"] == "phys" and r["qubit"] == q
                         and r["bit"] == bit and r["K"] == K]
                    vals.append(f"{v[0]['err']:>10.4f}" if v else f"{'—':>10s}")
                print(f"{K:>4d}", *vals)
    print(f"\n已写 {OUT / f'summary_{TAG}.json'}")


if __name__ == "__main__":
    main()
