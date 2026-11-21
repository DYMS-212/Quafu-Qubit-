"""counts 分析 + 统计严谨性工具。

count 的 bit 串序 = c_{n-1}..c1 c0(最右 = c0)。
本项目铁律:headline 数字必须带置信区间(Wilson 单比例 / bootstrap 通用)。
"""
from __future__ import annotations

import math
import random


def merge_counts(counts_list) -> dict:
    out: dict = {}
    for c in counts_list:
        for k, v in (c or {}).items():
            out[k] = out.get(k, 0) + v
    return out


def shots_of(count: dict) -> int:
    return sum(count.values())


def bit_p1(count: dict, n: int) -> list[float]:
    """每个逻辑位 i(从右数)读到 1 的概率。"""
    ones = [0] * n
    tot = 0
    for k, v in count.items():
        tot += v
        ks = k.zfill(n)
        for i in range(n):
            if ks[n - 1 - i] == "1":
                ones[i] += v
    tot = tot or 1
    return [o / tot for o in ones]


def bits_of(key: str, n: int) -> list[int]:
    """bit 串 → [b0, b1, ..., b_{n-1}](按逻辑位序,b0 = 最右)。"""
    ks = key.zfill(n)
    return [int(ks[n - 1 - i]) for i in range(n)]


def assignment_from_p1(p1_prep0: list[float], p1_prep1: list[float]) -> list[dict]:
    """由 P(1|备0)、P(1|备1) 算每比特读出指标(F0/F1/eta/bias)。"""
    out = []
    for a, b in zip(p1_prep0, p1_prep1):
        F0, F1 = 1 - a, b
        out.append({"F0": F0, "F1": F1, "eta": F0 + F1 - 1,
                    "bias": F0 - F1, "F_mean": (F0 + F1) / 2})
    return out


# ------------------------------------------------------------- 置信区间
def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """二项比例 k/n 的 Wilson 95% CI。"""
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def bootstrap_counts(count: dict, rng: random.Random) -> dict:
    """按多项分布重采样一份 counts(shots 不变)。"""
    keys = list(count.keys())
    weights = [count[k] for k in keys]
    n = sum(weights)
    draws = rng.choices(keys, weights=weights, k=n)
    out: dict = {}
    for k in draws:
        out[k] = out.get(k, 0) + 1
    return out


def bootstrap_ci(count: dict, stat_fn, B: int = 300, seed: int = 0,
                 z_lo: float = 2.5, z_hi: float = 97.5) -> dict:
    """对 stat_fn(count)->float 做 bootstrap,返回 {value, lo, hi}。"""
    rng = random.Random(seed)
    val = stat_fn(count)
    samples = sorted(stat_fn(bootstrap_counts(count, rng)) for _ in range(B))
    lo = samples[int(len(samples) * z_lo / 100)]
    hi = samples[min(len(samples) - 1, int(len(samples) * z_hi / 100))]
    return {"value": val, "lo": lo, "hi": hi, "B": B}


# ------------------------------------------------------------- 跨轮/跨日聚合
def pool_ci(values, B: int = 2000, seed: int = 0) -> dict:
    """多轮/多日标量值列表 → 聚合统计:mean、bootstrap 95% CI、best(min)、
    worst(max)、median、sd、n。把 headline 的多轮重复汇成带误差棒的数字
    (单轮 Wilson 只反映 shot 噪声,跨轮 pool 才反映 job-to-job + 日间漂移)。
    注:best 取 min(错误率越小越好);若指标越大越好由调用方读 worst/自行取。"""
    vals = [v for v in values if v is not None]
    n = len(vals)
    if n == 0:
        return {"n": 0}
    m = sum(vals) / n
    sv = sorted(vals)
    med = sv[n // 2] if n % 2 else (sv[n // 2 - 1] + sv[n // 2]) / 2
    if n == 1:
        return {"mean": m, "sd": 0.0, "lo": m, "hi": m, "n": 1,
                "best": m, "worst": m, "median": m}
    sd = (sum((v - m) ** 2 for v in vals) / (n - 1)) ** 0.5
    rng = random.Random(seed)
    means = sorted(sum(rng.choice(vals) for _ in range(n)) / n for _ in range(B))
    return {"mean": m, "sd": sd,
            "lo": means[int(B * 0.025)], "hi": means[min(B - 1, int(B * 0.975))],
            "n": n, "best": sv[0], "worst": sv[-1], "median": med}


# ------------------------------------------------------------- 深度扫描的每层错误率
def fit_decay_per_layer(pairs) -> dict:
    """(k, p) 点列 → 每层错误率 r:拟合 ln(1-2p) = a + k·ln(1-2r)(RB 式衰减)。

    比"交叉点"稳健:交叉与否由两臂截距(制备/读出偏置)决定、随日漂移,
    斜率只反映每层门质量。p≥0.5 的点(已饱和)自动剔除;返回 r 及其 95% CI
    (斜率标准误传播)与截距 a(≈制备+读出损耗)。点数 <3 时返回 {"n": ...}。"""
    pts = [(k, p) for k, p in pairs if p is not None and 0 <= p < 0.5]
    if len(pts) < 3:
        return {"n": len(pts)}
    xs = [k for k, _ in pts]
    ys = [math.log(1 - 2 * p) for _, p in pts]
    n = len(pts)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return {"n": n}
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    a = my - b * mx
    resid = [y - (a + b * x) for x, y in zip(xs, ys)]
    se_b = math.sqrt((sum(e * e for e in resid) / max(1, n - 2)) / sxx)
    to_r = lambda s: (1 - math.exp(s)) / 2  # noqa: E731
    return {"r": to_r(b), "r_lo": to_r(b + 1.96 * se_b), "r_hi": to_r(b - 1.96 * se_b),
            "slope": b, "intercept": a, "prep_offset": to_r(a), "n": n}
