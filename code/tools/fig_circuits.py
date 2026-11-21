"""实验线路图集:qiskit mpl 渲染器直接画我们提交的 OpenQASM(零额度)。

画的就是真实提交线路(builder 产出的 QASM 原文,非示意重绘);**线路导轨以
物理比特号(Q67 等)标注**,取自各实验最近一轮 summary 的实际选址(每日体检后
自动更新);拿不到映射时回退逻辑编号。产出 figures/<日期>/circuit_*.png/svg/pdf。

依赖 vendor_viz/(qiskit+pylatexenc,Mac 本地,pip install --target vendor_viz)。
运行(Mac): ./qrun tools/fig_circuits.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))
sys.path.append(str(PROJ / "vendor_viz"))     # qiskit 专用,排最后不干扰主依赖

import matplotlib.pyplot as plt  # noqa: E402
from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, qasm2  # noqa: E402

from qfab import repcode  # noqa: E402
from tools import fig_style as st  # noqa: E402


def _load(modname: str):
    spec = importlib.util.spec_from_file_location(
        modname, PROJ / "experiments" / f"{modname}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def latest(exp: str):
    """最近一轮 summary(取物理选址用);无则 None。"""
    files = sorted((PROJ / "results" / exp).glob("summary_2*.json"))
    return json.loads(files[-1].read_text()) if files else None


def _get(d, *keys, default=None):
    for k in keys:
        if not isinstance(d, dict):
            return default
        d = d.get(k)
    return d if d is not None else default


def relabel(qc: QuantumCircuit, phys: list[int]) -> QuantumCircuit:
    """把逻辑导轨 q0..qn 重挂到以物理比特号命名的寄存器(Q67 …)。
    顺序约定 = target_qubits:第 i 条导轨的物理比特 = phys[i](全套 builder 一致)。"""
    regs = [QuantumRegister(1, name=f"Q{q}") for q in phys]
    new = QuantumCircuit(*regs, ClassicalRegister(qc.num_clbits, "c"))
    new.compose(qc, qubits=range(qc.num_qubits),
                clbits=range(qc.num_clbits), inplace=True)
    return new


CX_NOTE = "cx=OpenQASM 通用写法,平台编译为原生门集 h/rx/ry/rz/cz(执行视图另见 *_native 图)"


def draw(name: str, qasm: str, title: str, note: str = "",
         phys: list[int] | None = None, scale: float = 0.9):
    qc = qasm2.loads(qasm)
    if phys and len(phys) == qc.num_qubits:
        qc = relabel(qc, phys)
        note = (note + ";" if note else "") + "比特号=最近一轮实际选址(每日自动更新)"
    fig = qc.draw("mpl", scale=scale, fold=26)
    fig.suptitle(title, fontsize=13, color=st.INK, y=1.02)
    if "cx" in qasm:
        note = (note + "。" if note else "") + CX_NOTE
    if note:
        fig.text(0.01, -0.02, note, fontsize=8.5, color=st.MUTED)
    st.save(fig, f"circuit_{name}", formats=("png", "svg", "pdf"))


# ------------------------------------------------------------------ 执行视图
def latest_transpiled(exp: str, label_sub: str, exclude: str = "\x00") -> str | None:
    """最近一轮记录里该 label 的平台返回编译线路(真实执行视图)。"""
    root = PROJ / "results" / exp
    if not root.exists():
        return None
    for d in sorted((x for x in root.iterdir() if x.is_dir()), reverse=True):
        for f in sorted(d.glob("*.json")):
            try:
                r = json.loads(f.read_text())
            except Exception:  # noqa: BLE001
                continue
            lab = r.get("label", "")
            if label_sub in lab and exclude not in lab:
                t = (r.get("result") or {}).get("transpiled")
                if t:
                    return t
    return None


def compact(t: str) -> tuple[str, list[int]]:
    """把 156 线的 transpiled 压缩到只含用到的比特(保持顺序),返回 (qasm, 物理序)。
    "用到" = 出现在门/测量行;仅出现在 barrier 里的比特(空导轨)剔除。"""
    import re
    used = sorted({int(m) for ln in t.splitlines()
                   if not ln.strip().startswith("barrier")
                   for m in re.findall(r"q\[(\d+)\]", ln)})
    idx = {q: i for i, q in enumerate(used)}
    lines = []
    for ln in t.splitlines():
        if ln.startswith("qreg"):
            lines.append(f"qreg q[{len(used)}];")
        elif ln.strip().startswith("barrier"):
            qs = [int(m) for m in re.findall(r"q\[(\d+)\]", ln)]
            keep = [q for q in qs if q in idx]
            if not qs:
                lines.append(ln)              # barrier q;
            elif keep:
                lines.append("barrier " + ",".join(f"q[{idx[q]}]" for q in keep) + ";")
        else:
            lines.append(re.sub(r"q\[(\d+)\]",
                                lambda mm: f"q[{idx[int(mm.group(1))]}]", ln))
    return "\n".join(lines), used


def draw_native(name: str, exp: str, label_sub: str, title: str,
                exclude: str = "\x00", scale: float = 0.85):
    t = latest_transpiled(exp, label_sub, exclude)
    if not t:
        print(f"跳过 {name}_native(暂无 transpiled 记录)")
        return
    qasm, used = compact(t)
    qc = relabel(qasm2.loads(qasm), used)
    fig = qc.draw("mpl", scale=scale, fold=30)
    fig.suptitle(title + " · 执行视图(平台返回的编译线路)",
                 fontsize=13, color=st.INK, y=1.02)
    fig.text(0.01, -0.02, "所见即所执行:平台编译器输出,已映射物理比特;"
                          "原生门集 h/rx/ry/rz/cz(个别 cx 为平台编译层保留,底层再分解)",
             fontsize=8.5, color=st.MUTED)
    st.save(fig, f"circuit_{name}_native", formats=("png", "svg", "pdf"))


def main():
    st.apply_style()
    e10 = _load("e10_code422")
    e11 = _load("e11_encoded_memory")
    e18 = _load("e18_t2_swap")
    e19 = _load("e19_code422_gates")

    # 最近一轮的真实物理选址(缺则回退逻辑编号)
    s05 = latest("e05_repcode")
    rows05 = s05 if isinstance(s05, list) else _get(s05, "rows", default=[])
    chain3 = next((r["chain"] for r in rows05
                   if r.get("code") == "bit" and r.get("d") == 3), None)
    q10 = _get(latest("e10_code422"), "provenance", "qubits", default={})
    star = (q10.get("stars") or [None])[0]
    fusion = (q10.get("fusion_paths") or [None])[0]
    q11 = _get(latest("e11_memory"), "provenance", "qubits", default={})
    mem = ((q11.get("chains") or {}).get("3", []) + (q11.get("load_edge") or [])) or None
    star19 = _get(latest("e19_gates"), "provenance", "qubits", "star")
    lay18 = _get(latest("e18_t2swap"), "layout")

    draw("repcode_d3",
         repcode.build_t1(3, logical_bit=1),
         "重复码(码距 3):编码 → 稳定子提取 → 末端测量",
         "导轨 1/3/5 = 数据比特,2/4 = 检验比特;逻辑 |1> = |111>",
         phys=chain3)
    draw("repcode_d3_phase",
         repcode.build_t1(3, logical_bit=1, phase=True),
         "相位型重复码(码距 3):H 夹层让退相位成为可纠错误",
         "数据比特在 X 基存储;提取前转回 Z 基", phys=chain3)
    draw("code422_star",
         e10.circ_422_z(0, 0),
         "[[4,2,2]] 星形制备:GHZ 编码 |00>_L + 横向读出",
         "度 3 星形放置(首导轨=星心);末端后选稳定子奇偶即检错", phys=star)
    draw("code422_fusion",
         e10.circ_fusion_z(0, 0),
         "[[4,2,2]] 验证式制备:两个 Bell 对 + 奇偶校验比特融合",
         "第 3 导轨 = 验证比特(后选 0 = 投影成逻辑态 + 制备检错),其余为码比特",
         phys=fusion)
    draw("memory_d3K4",
         e11.build_encoded(3, 1, 4)[0],
         "编码存储:编码 → K 层负载(存储期,错误真实发生)→ 提取",
         "末两条导轨 = 远端负载边(逐层 barrier 同步计时);图示 K=4,实验扫 K=0..32",
         phys=mem)
    draw("memory_d3K4_phase",
         e11.build_encoded(3, 1, 4, phase=True)[0],
         "相位编码存储:存 |−>_L,存储期积累的是 T2 退相位",
         "同一负载与提取;对'重复码=经典冗余'质疑的正面回答", phys=mem)
    draw("gates_k2",
         e19.circ_encoded(2),
         "逻辑门序列:横向 H × k(图示 k=2)",
         "[[4,2,2]] 的横向 H = 双逻辑 H·SWAP;层间 barrier 防编译器对消",
         phys=star19)
    draw("t2swap",
         e18.circ("none"),
         "两轮稳定子测量:SWAP 换出合成(无中途测量)",
         "第 3 导轨=检验比特,round1 结果经 SWAP 存入第 4 导轨,原位复用测 round2",
         phys=lay18)

    # 执行视图(真实 transpiled,原生门+物理映射;评委问"机器没有 CNOT"的标准答案)
    draw_native("code422_star", "e10_code422", "_Z_00", exclude="v",
                title="[[4,2,2]] 星形制备")
    draw_native("code422_fusion", "e10_code422", "_Z_00", exclude="c",
                title="[[4,2,2]] 验证式制备")
    draw_native("memory_d3_phase", "e11_memory", "d3P_L0_K8",
                title="相位编码存储(d=3, K=8)")
    draw_native("logical_bell", "e24_logical_bell", "logical_r0",
                title="逻辑 Bell 态制备(transversal S)")
    print("线路图集完成 →", st.outdir())


if __name__ == "__main__":
    main()
