"""repetition 码(单轮 T=1,fresh ancilla,CZ 原生机可跑)。

线性链布局 chain = [D0, A0, D1, A1, ..., D_{d-1}]:
  data  D_k -> 逻辑位 2k;ancilla A_k -> 逻辑位 2k+1,检验 Z_{D_k}Z_{D_{k+1}}
bit-flip 码:|0>_L=|0..0>,|1>_L=|1..1>;phase-flip 码:同构,前后加 H 夹层。
"""
from __future__ import annotations

from .qasm import QASM


def chain_layout(data_phys: list[int], anc_phys: list[int]) -> list[int]:
    d = len(data_phys)
    assert len(anc_phys) == d - 1
    chain = []
    for k in range(d):
        chain.append(data_phys[k])
        if k < d - 1:
            chain.append(anc_phys[k])
    return chain


def build_t1(d: int, logical_bit: int = 0, inject_data: list[int] | None = None,
             phase: bool = False, idle_layers: int = 0) -> str:
    """单轮 repetition 码线路(OpenQASM 2.0),配合 target_qubits=chain_layout(...)。

    phase=True: phase-flip 码(data 上 H 夹层,注入错误为 Z)。
    idle_layers>0: 编码后、syndrome 提取前插入 K 层 barrier 空转(编码存储实验用,
                   真 idle 需另配负载比特,见实验脚本)。
    """
    inject = inject_data or []
    n = 2 * d - 1
    qc = QASM(n)
    data = [2 * k for k in range(d)]
    # 制备逻辑态
    if logical_bit:
        for q in data:
            qc.x(q)
    if phase:
        for q in data:
            qc.h(q)
    qc.barrier()
    # idle 存储期(错误在这里真实发生)
    for _ in range(idle_layers):
        qc.barrier()
    # 注入错误
    for k in inject:
        if phase:
            qc.z(2 * k)
        else:
            qc.x(2 * k)
    qc.barrier()
    # phase 码:检验前把 data 转回 Z 基(H),等价于在 X 基做 ZZ 检验
    if phase:
        for q in data:
            qc.h(q)
    # syndrome 提取:A_k = D_k ⊕ D_{k+1}
    for k in range(d - 1):
        qc.cx(2 * k, 2 * k + 1)
        qc.cx(2 * (k + 1), 2 * k + 1)
    qc.barrier()
    qc.measure_all()
    return qc.qasm()


# ------------------------------------------------------------- 基础解码
def split_bits(key: str, d: int) -> tuple[list[int], list[int]]:
    """bit 串 → (data[d], syndrome[d-1])。逻辑位 2k=data,2k+1=ancilla。"""
    n = 2 * d - 1
    ks = key.zfill(n)
    data = [int(ks[n - 1 - 2 * k]) for k in range(d)]
    synd = [int(ks[n - 1 - (2 * k + 1)]) for k in range(d - 1)]
    return data, synd


def majority(bits: list[int]) -> int:
    return 1 if 2 * sum(bits) > len(bits) else 0


def logical_error_majority(count: dict, d: int, prepared_bit: int) -> float:
    """多数表决逻辑错误率(忽略 syndrome 位)。"""
    tot = err = 0
    for k, v in count.items():
        data, _ = split_bits(k, d)
        if majority(data) != prepared_bit:
            err += v
        tot += v
    return err / (tot or 1)
