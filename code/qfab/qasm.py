"""OpenQASM 2.0 构造器(平台编译器会把 cx/h 等转译到原生门 h/rx/ry/rz/cz)。

特性:
  - pauli(p, q) 便捷门,twirl_cz() 对单个 CZ 做 Pauli twirling 包裹
  - idle_block() 用 barrier 分层,配合负载比特造真 idle(编码存储实验用)
"""
from __future__ import annotations

import math
import random

PI = math.pi

# CZ twirling 查找表:CZ·(P_c⊗P_t) = (P'_c⊗P'_t)·CZ(相位差全局,可忽略)
# key = 输入两侧 Pauli,value = 输出两侧 Pauli
_CZ_TWIRL = {
    ("I", "I"): ("I", "I"), ("I", "X"): ("Z", "X"), ("I", "Y"): ("Z", "Y"), ("I", "Z"): ("I", "Z"),
    ("X", "I"): ("X", "Z"), ("X", "X"): ("Y", "Y"), ("X", "Y"): ("Y", "X"), ("X", "Z"): ("X", "I"),
    ("Y", "I"): ("Y", "Z"), ("Y", "X"): ("X", "Y"), ("Y", "Y"): ("X", "X"), ("Y", "Z"): ("Y", "I"),
    ("Z", "I"): ("Z", "I"), ("Z", "X"): ("I", "X"), ("Z", "Y"): ("I", "Y"), ("Z", "Z"): ("Z", "Z"),
}


class QASM:
    def __init__(self, nq: int, nc: int | None = None):
        self.nq = nq
        self.nc = nc if nc is not None else nq
        self.ops: list[str] = []

    # ---- 单比特门 ----
    def h(self, q):
        self.ops.append(f"h q[{q}];"); return self

    def x(self, q):
        self.ops.append(f"x q[{q}];"); return self

    def y(self, q):
        self.ops.append(f"y q[{q}];"); return self

    def z(self, q):
        self.ops.append(f"z q[{q}];"); return self

    def s(self, q):
        self.ops.append(f"s q[{q}];"); return self

    def sdg(self, q):
        self.ops.append(f"sdg q[{q}];"); return self

    def rx(self, theta, q):
        self.ops.append(f"rx({theta}) q[{q}];"); return self

    def ry(self, theta, q):
        self.ops.append(f"ry({theta}) q[{q}];"); return self

    def rz(self, theta, q):
        self.ops.append(f"rz({theta}) q[{q}];"); return self

    def pauli(self, p: str, q):
        if p in ("X", "Y", "Z"):
            self.ops.append(f"{p.lower()} q[{q}];")
        return self  # "I" 不发指令

    # ---- 两比特门 ----
    def cx(self, c, t):
        self.ops.append(f"cx q[{c}],q[{t}];"); return self

    def cz(self, c, t):
        self.ops.append(f"cz q[{c}],q[{t}];"); return self

    def twirl_cz(self, c, t, rng: random.Random):
        """Pauli-twirled CZ:随机 P 包裹,逻辑上仍等价于 CZ。"""
        pc, pt = rng.choice(list(_CZ_TWIRL))
        qc, qt = _CZ_TWIRL[(pc, pt)]
        self.pauli(pc, c).pauli(pt, t)
        self.cz(c, t)
        self.pauli(qc, c).pauli(qt, t)
        return self

    # ---- 结构 ----
    def barrier(self, *qs):
        if qs:
            self.ops.append("barrier " + ",".join(f"q[{q}]" for q in qs) + ";")
        else:
            self.ops.append("barrier q;")
        return self

    def measure(self, q, c):
        self.ops.append(f"measure q[{q}] -> c[{c}];"); return self

    def measure_all(self):
        for i in range(min(self.nq, self.nc)):
            self.measure(i, i)
        return self

    # ---- 输出 ----
    def qasm(self) -> str:
        head = [
            "OPENQASM 2.0;",
            'include "qelib1.inc";',
            f"qreg q[{self.nq}];",
            f"creg c[{self.nc}];",
        ]
        return "\n".join(head + self.ops) + "\n"

    def __str__(self):
        return self.qasm()


def bell_pair(basis: str = "ZZ") -> QASM:
    """两比特 Bell 态 (|00>+|11>)/√2,basis='ZZ'|'XX' 决定测量基。"""
    qc = QASM(2)
    qc.h(0).cx(0, 1)
    if basis == "XX":
        qc.h(0).h(1)
    qc.barrier().measure_all()
    return qc
