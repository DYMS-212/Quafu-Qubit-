"""读出误差缓解:张量化 assignment 矩阵反演(标准方法,自测矩阵驱动)。

模型:测量分布 p_meas = (⊗_i A_i) p_true,A_i = [[F0, 1-F1], [1-F0, F1]](列=真值,行=读出)。
缓解:p_true ≈ (⊗_i A_i^{-1}) p_meas,逐比特作用(O(n·2^n),n≤15 无压力),
     负值截断后重归一(标注截断量,诚实报告)。
"""
from __future__ import annotations


def correct_counts(count: dict, fs: list[tuple[float, float]]) -> tuple[dict, float]:
    """count: bitstring→次数(位序 c_{n-1}..c0,逻辑位 i 从右数);
    fs[i] = (F0, F1) 为逻辑位 i 的读出保真(按 target_qubits 映射好再传入)。
    返回 (缓解后的伪计数 dict, 负值截断占比)。"""
    n = len(fs)
    dim = 1 << n
    v = [0.0] * dim
    tot = 0
    for k, c in count.items():
        v[int(k.zfill(n), 2)] = float(c)
        tot += c
    if not tot:
        return {}, 0.0
    for i in range(n):  # 逻辑位 i(从右数)= 整数位 i
        F0, F1 = fs[i]
        det = F0 + F1 - 1
        if det <= 0.05:   # 读出接近随机,反演病态 → 跳过该位
            continue
        a, b = F1 / det, (F1 - 1) / det   # A^{-1} = 1/det * [[F1, F1-1],[F0-1, F0]]
        c2, d2 = (F0 - 1) / det, F0 / det
        step = 1 << i
        for base in range(0, dim, step << 1):
            for off in range(base, base + step):
                v0, v1 = v[off], v[off + step]
                v[off] = a * v0 + b * v1
                v[off + step] = c2 * v0 + d2 * v1
    clipped = sum(-x for x in v if x < 0)
    v = [max(x, 0.0) for x in v]
    s = sum(v) or 1.0
    scale = tot / s
    out = {format(idx, f"0{n}b"): x * scale for idx, x in enumerate(v) if x > 1e-9}
    return out, clipped / tot


def fs_for_qubits(qubits: list[int], cal: dict[int, dict],
                  default: tuple[float, float] = (0.97, 0.93)) -> list[tuple[float, float]]:
    """按 target_qubits 顺序取每位 (F0,F1);cal = {物理比特: {F0,F1,...}}(e02 产出)。"""
    return [(cal[q]["F0"], cal[q]["F1"]) if q in cal else default for q in qubits]
