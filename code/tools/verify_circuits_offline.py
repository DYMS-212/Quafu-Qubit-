"""离线线路验证器:statevector 精确模拟实验脚本产出的 QASM,验证逻辑正确性。

零额度。解析的是 build/circ 函数**真实输出**的 QASM 串(非重新实现),
过一遍与真机分析完全相同的 eval 函数,防"线路对/分析错"或反之。
覆盖:e10 星形(回归)、e10 融合制备(新)、e11 bit/phase 编码存储(phase 新)。

运行: ./qrun tools/verify_circuits_offline.py
"""
from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

import numpy as np

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))

from qfab import repcode  # noqa: E402

_GATE1 = {
    "h": np.array([[1, 1], [1, -1]]) / np.sqrt(2),
    "x": np.array([[0, 1], [1, 0]], dtype=complex),
    "y": np.array([[0, -1j], [1j, 0]]),
    "z": np.array([[1, 0], [0, -1]], dtype=complex),
    "s": np.array([[1, 0], [0, 1j]]),
    "sdg": np.array([[1, 0], [0, -1j]]),
}


def simulate(qasm: str) -> dict:
    """QASM → 理想测量分布 {c串: 概率}(c 串右端=c0,与平台 count 位序一致)。"""
    nq = nc = 0
    psi = None
    cmap: dict[int, int] = {}          # c[j] = q[i]
    for line in qasm.splitlines():
        line = line.strip()
        if not line or line.startswith(("OPENQASM", "include", "barrier", "//")):
            continue
        m = re.match(r"qreg q\[(\d+)\];", line)
        if m:
            nq = int(m.group(1))
            psi = np.zeros(2 ** nq, dtype=complex)
            psi[0] = 1.0
            continue
        m = re.match(r"creg c\[(\d+)\];", line)
        if m:
            nc = int(m.group(1))
            continue
        m = re.match(r"(h|x|y|z|s|sdg) q\[(\d+)\];", line)
        if m:
            g, q = m.group(1), int(m.group(2))
            v = psi.reshape(2 ** (nq - q - 1), 2, 2 ** q)
            psi = np.einsum("ab,ibj->iaj", _GATE1[g], v).reshape(-1)
            continue
        m = re.match(r"(rx|ry|rz)\(([-+0-9.eE]+)\) q\[(\d+)\];", line)
        if m:
            g, th, q = m.group(1), float(m.group(2)), int(m.group(3))
            c_, s_ = np.cos(th / 2), np.sin(th / 2)
            u = {"rx": np.array([[c_, -1j * s_], [-1j * s_, c_]]),
                 "ry": np.array([[c_, -s_], [s_, c_]], dtype=complex),
                 "rz": np.array([[np.exp(-0.5j * th), 0],
                                 [0, np.exp(0.5j * th)]])}[g]
            v = psi.reshape(2 ** (nq - q - 1), 2, 2 ** q)
            psi = np.einsum("ab,ibj->iaj", u, v).reshape(-1)
            continue
        m = re.match(r"(cx|cz) q\[(\d+)\],q\[(\d+)\];", line)
        if m:
            g, c, t = m.group(1), int(m.group(2)), int(m.group(3))
            idx = np.arange(2 ** nq)
            both = ((idx >> c) & 1).astype(bool)
            if g == "cz":
                sign = np.where(both & (((idx >> t) & 1) == 1), -1.0, 1.0)
                psi = psi * sign
            else:
                new = psi.copy()
                sel = idx[both]
                new[sel] = psi[sel ^ (1 << t)]
                psi = new
            continue
        m = re.match(r"measure q\[(\d+)\] -> c\[(\d+)\];", line)
        if m:
            cmap[int(m.group(2))] = int(m.group(1))
            continue
        raise ValueError(f"验证器不识别的 QASM 行: {line}")
    probs = np.abs(psi) ** 2
    out: dict[str, float] = {}
    for i in np.nonzero(probs > 1e-12)[0]:
        cval = 0
        for cj, qi in cmap.items():
            cval |= ((int(i) >> qi) & 1) << cj
        key = format(cval, f"0{nc}b")
        out[key] = out.get(key, 0.0) + float(probs[i])
    return out


def _load(modname: str):
    spec = importlib.util.spec_from_file_location(
        modname, PROJ / "experiments" / f"{modname}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def check(name: str, ok: bool, detail: str = ""):
    print(f"  {'✅' if ok else '❌'} {name}" + (f"  {detail}" if detail else ""))
    return ok


def main():
    all_ok = True
    e10 = _load("e10_code422")

    print("== e10 星形 [[4,2,2]](回归)==")
    for basis, builder in (("Z", e10.circ_422_z), ("X", e10.circ_422_x)):
        for m in (0, 1):
            for n in (0, 1):
                ev = e10.eval_422(simulate(builder(m, n)), basis, m, n)
                all_ok &= check(f"star {basis}{m}{n}",
                                abs(ev["acceptance"] - 1) < 1e-9
                                and ev["ps_err_perq"] < 1e-9
                                and ev["raw_err_perq"] < 1e-9)

    print("== e10 融合制备(新)==")
    for basis, builder in (("Z", e10.circ_fusion_z), ("X", e10.circ_fusion_x)):
        for m in (0, 1):
            for n in (0, 1):
                ev = e10.eval_fusion(simulate(builder(m, n)), basis, m, n)
                all_ok &= check(
                    f"fusion {basis}{m}{n}",
                    abs(ev["v_acceptance"] - 0.5) < 1e-9
                    and abs(ev["acceptance"] - 0.5) < 1e-9
                    and ev["ps_err_perq"] < 1e-9 and ev["raw_err_perq"] < 1e-9,
                    f"v_acc={ev['v_acceptance']:.3f}")

    # 融合容错账本:逐位置注入单 X 故障。码比特/Bell/CX 故障必须全部被
    # v 或 S_Z 剔除(err=0);**唯一允许的泄漏 = v 自身在首个 CX 之后的翻转**
    # (选错融合分支 ≡ X̄₁,err=0.5)——这是方案的已知残余暴露,写进 e10 docstring。
    print("== 融合容错账本(注入单 X 故障;v 窗残余=预期泄漏)==")
    base = e10.circ_fusion_z(0, 0)
    lines = base.splitlines()
    v_window = False                       # v 的首个 CX 之后 → 分支比特已在 v 上
    for i, ln in enumerate(lines):
        if not ln.startswith(("cx", "h ")):
            continue
        if ln.startswith("cx") and ",q[2]" in ln:
            v_window = True
        for q in sorted({int(x) for x in re.findall(r"q\[(\d+)\]", ln)}):
            fv = lines[:i + 1] + [f"x q[{q}];"] + lines[i + 1:]
            ev = e10.eval_fusion(simulate("\n".join(fv)), "Z", 0, 0)
            expected_leak = (q == 2 and v_window)
            if expected_leak:
                ok = abs(ev["ps_err_perq"] - 0.5) < 1e-9
                tag = "预期泄漏(v 窗)"
            else:
                ok = ev["ps_err_perq"] < 1e-9
                tag = "必须被检"
            all_ok &= check(f"故障[{ln.strip()} 后 X q{q}] {tag}", ok,
                            f"acc={ev['acceptance']:.3f} err={ev['ps_err_perq']:.4f}")

    print("== e11 编码存储 bit/phase(phase 新)==")
    e11 = _load("e11_encoded_memory")
    for phase in (False, True):
        for d in (3, 5):
            for bit in (0, 1):
                for K in (0, 2):
                    qasm, _ = e11.build_encoded(d, bit, K, phase=phase)
                    cnt = simulate(qasm)
                    n = 2 * d - 1
                    sub: dict = {}
                    for key, v in cnt.items():
                        kk = key.zfill(n + 2)[2:]
                        sub[kk] = sub.get(kk, 0) + v
                    err = repcode.logical_error_majority(sub, d, bit)
                    all_ok &= check(
                        f"enc {'phase' if phase else 'bit'} d{d} L{bit} K{K}",
                        err < 1e-9)
    for phase in (False, True):
        for bit in (0, 1):
            for K in (0, 2):
                qasm, _ = e11.build_phys(bit, K, 0, phase=phase)
                cnt = simulate(qasm)
                p1 = sum(v for k, v in cnt.items() if k[-1] == "1")
                all_ok &= check(
                    f"phys {'phase' if phase else 'bit'} L{bit} K{K}",
                    abs(p1 - bit) < 1e-9)

    print("\n" + ("全部通过 ✅" if all_ok else "有失败 ❌,勿上真机"))
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
