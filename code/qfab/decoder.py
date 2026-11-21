"""校准感知 MWPM 解码(pymatching),针对单轮 T=1 repetition + 末端测量。

物理模型(重要,报告要讲):
  设 prep 逻辑态 L,物理数据错 e_k(提取前)、数据读出翻转 m^d_k、ancilla 读出翻转 m^a_k:
    观测 D_k = L ⊕ e_k ⊕ m^d_k
    观测 s_k = e_k ⊕ e_{k+1} ⊕ m^a_k
  detector: det_k = s_k ⊕ D_k ⊕ D_{k+1} = m^a_k ⊕ m^d_k ⊕ m^d_{k+1}   ← e 完全抵消!
  ⇒ 在此设定下 syndrome 的信息量 = 对"读出/测量类错误"的一致性校验;
    物理数据错 e 由(修正后的)多数表决处理。
  ⇒ 解码 = 在 detector 图上 MWPM 判别 m^d vs m^a,把判为 m^d 的位翻回来再表决。
    权重来自自测 assignment matrix ⇒ "校准感知解码",自校准主线第三处闭环。

匹配图(链 + 边界):
  节点 = det_0..det_{d-2} + 边界
  数据位 j 的读出错: 连 det_{j-1}—det_j(端点连边界),fault_id=j(要翻数据位)
  ancilla k 的读出错: det_k—边界,无 fault(不影响逻辑)
"""
from __future__ import annotations

import math

from .repcode import split_bits, majority


def _w(p: float) -> float:
    p = min(max(p, 1e-6), 0.499999)
    return math.log((1 - p) / p)


def build_matching(d: int, p_data: list[float], p_anc: list[float]):
    """p_data[d]、p_anc[d-1] = 各位读出翻转概率(uniform 或自测)。"""
    import pymatching
    assert len(p_data) == d and len(p_anc) == d - 1
    m = pymatching.Matching()
    boundary = d - 1  # detector 节点 0..d-2,边界节点 d-1
    for j in range(d):  # 数据读出错
        u = j - 1 if j - 1 >= 0 else boundary
        v = j if j <= d - 2 else boundary
        m.add_edge(u, v, fault_ids={j}, weight=_w(p_data[j]),
                   merge_strategy="smallest-weight")
    for k in range(d - 1):  # ancilla 读出错
        m.add_edge(k, boundary, fault_ids=set(), weight=_w(p_anc[k]),
                   merge_strategy="smallest-weight")
    m.set_boundary_nodes({boundary})
    return m


def decode_shot(matching, data: list[int], synd: list[int]) -> list[int]:
    """返回修正后的 data 位。"""
    d = len(data)
    det = [synd[k] ^ data[k] ^ data[k + 1] for k in range(d - 1)]
    if not any(det):
        return data
    corr = matching.decode(det)  # 长度 d 的 0/1 数组(fault_ids=数据位)
    return [b ^ int(c) for b, c in zip(data, corr)]


def logical_error(count: dict, d: int, prepared_bit: int,
                  p_data: list[float], p_anc: list[float],
                  soft: bool = False) -> float:
    """MWPM 修正读出错误后表决的逻辑错误率。
    soft=True: 修正后用校准加权表决(权重 ln((1-p)/p))而非简单多数。"""
    matching = build_matching(d, p_data, p_anc)
    wts = [_w(p) for p in p_data]
    cache: dict = {}
    tot = err = 0
    for key, v in count.items():
        if key not in cache:
            data, synd = split_bits(key, d)
            fixed = decode_shot(matching, data, synd)
            if soft:
                score = sum(w * (1 if b else -1) for w, b in zip(wts, fixed))
                cache[key] = 1 if score > 0 else 0
            else:
                cache[key] = majority(fixed)
        if cache[key] != prepared_bit:
            err += v
        tot += v
    return err / (tot or 1)


def detector_rates(count: dict, d: int) -> list[float]:
    """经验 detector 触发率 P(det_k=1)。det_k ≈ p_anc_k + p_data_k + p_data_{k+1}(小 p 近似),
    可用于交叉拟合 ancilla 有效错误率(含 CX 门错,readout-only 先验严重低估它)。"""
    tot = 0
    fired = [0] * (d - 1)
    for key, v in count.items():
        data, synd = split_bits(key, d)
        tot += v
        for k in range(d - 1):
            if synd[k] ^ data[k] ^ data[k + 1]:
                fired[k] += v
    tot = tot or 1
    return [f / tot for f in fired]


def p_anc_from_detectors(fit_count: dict, d: int, p_data_fit: list[float],
                         floor: float = 5e-3) -> list[float]:
    """用另一份数据(如 L0 拟合、L1 解码)的 detector 率反解 ancilla 有效错误率。"""
    rates = detector_rates(fit_count, d)
    return [max(r - p_data_fit[k] - p_data_fit[k + 1], floor)
            for k, r in enumerate(rates)]


def _likelihood(data: list[int], synd: list[int], L: int,
                q: list[float], r: list[float], a: list[float]) -> float:
    """P(D,s | 逻辑值 L),对物理错 e 向量精确求和(链上 DP,O(d))。
    模型:D_k = L⊕e_k⊕m^d_k,s_k = e_k⊕e_{k+1}⊕m^a_k;
         P(e_k)=q_k,P(m^d_k)=r_k,P(m^a_k)=a_k,全独立。"""
    d = len(data)

    def p_obs(k, e):
        return r[k] if data[k] != (L ^ e) else 1 - r[k]

    alpha = [ (1 - q[0]) * p_obs(0, 0), q[0] * p_obs(0, 1) ]
    for k in range(1, d):
        nxt = [0.0, 0.0]
        for e2 in (0, 1):
            pe = q[k] if e2 else 1 - q[k]
            tot = 0.0
            for e1 in (0, 1):
                ps = a[k - 1] if synd[k - 1] != (e1 ^ e2) else 1 - a[k - 1]
                tot += alpha[e1] * ps
            nxt[e2] = tot * pe * p_obs(k, e2)
        alpha = nxt
    return alpha[0] + alpha[1]


def logical_error_map(count: dict, d: int, prepared_bit: int,
                      q: list[float], r: list[float], a: list[float]) -> float:
    """精确边缘化 MAP 解码:argmax_L P(D,s|L)。模型正确时可证明不劣于多数表决。
    q=物理翻转、r=数据读出翻转、a=ancilla 有效(读出+CX)——全部可自测/交叉拟合。"""
    cache: dict = {}
    tot = err = 0
    for key, v in count.items():
        if key not in cache:
            data, synd = split_bits(key, d)
            l1 = _likelihood(data, synd, 1, q, r, a)
            l0 = _likelihood(data, synd, 0, q, r, a)
            cache[key] = 1 if l1 > l0 else 0
        if cache[key] != prepared_bit:
            err += v
        tot += v
    return err / (tot or 1)


def logical_error_soft_majority(count: dict, d: int, prepared_bit: int,
                                p_data: list[float]) -> float:
    """不用 syndrome,仅校准加权表决(对照组:量化"软表决"单独的贡献)。"""
    wts = [_w(p) for p in p_data]
    tot = err = 0
    for key, v in count.items():
        data, _ = split_bits(key, d)
        score = sum(w * (1 if b else -1) for w, b in zip(wts, data))
        bit = 1 if score > 0 else 0
        if bit != prepared_bit:
            err += v
        tot += v
    return err / (tot or 1)
