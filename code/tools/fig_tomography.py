"""密度矩阵城市图:e23 层析结果 → qiskit plot_state_city 官方渲染(零额度)。

渲染内核 = qiskit.visualization.plot_state_city(用户定:负值的
"穿透 z=0 半透明地板下挂"处理最自然,单色即可,坐标轴已足够清楚)。
我们只做外壳:中文标题(机器/态/物理比特)、指标行(F/纯度/负性)、
低保真警示行、每图配套 CSV、缓解版、三机对比页。

产出 figures/<日期>/tomo_*.png/svg/pdf + rho CSV
运行(Mac,依赖 vendor_viz 的 qiskit): ./qrun tools/fig_tomography.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))
sys.path.append(str(PROJ / "vendor_viz"))

import matplotlib.pyplot as plt  # noqa: E402
from qiskit.quantum_info import DensityMatrix  # noqa: E402
from qiskit.visualization import plot_state_city  # noqa: E402

from tools import fig_style as st  # noqa: E402

R = PROJ / "results" / "e23_tomography"
CHIPS = ("Baihua", "Dongling", "Shenglian")


def latest(chip: str) -> dict | None:
    files = sorted(R.glob(f"summary_{chip}_*.json"))
    return json.loads(files[-1].read_text()) if files else None


def head(m: dict) -> str:
    s = f"F={m['fidelity']:.3f} · 纯度 {m['purity']:.3f}"
    if "negativity" in m:
        s += f" · 负性 {m['negativity']:.3f}"
    return s


def rho_of(m: dict) -> DensityMatrix:
    return DensityMatrix(np.array(m["rho_re"]) + 1j * np.array(m["rho_im"]))


def draw_state(chip: str, name: str, m: dict, suffix: str = ""):
    fig = plt.figure(figsize=(9.8, 4.8))
    axr = fig.add_subplot(1, 2, 1, projection="3d")
    axi = fig.add_subplot(1, 2, 2, projection="3d")
    plot_state_city(rho_of(m), ax_real=axr, ax_imag=axi)
    axr.set_title("Re ρ", fontsize=11, color=st.INK)
    axi.set_title("Im ρ", fontsize=11, color=st.INK)
    fig.subplots_adjust(left=0.0, right=1.0, wspace=0.04, top=0.84, bottom=0.08)
    tag = "(读出缓解后)" if suffix else ""
    qb = m.get("qubits") or []
    where = f"Q{'-Q'.join(str(q) for q in qb)}" if qb else "选址未记录"
    state_cn = "Bell 态" if name == "bell" else "GHZ 三比特态"
    fig.suptitle(f"{state_cn}断层成像{tag} · {chip} · {where}",
                 fontsize=12.5, y=0.98)
    fig.text(0.5, 0.885, head(m), ha="center", fontsize=9.5, color=st.INK2)
    if m["fidelity"] < 0.5:                    # 纠缠判据线:低于它=噪声支配态
        note = (f"[警示] F={m['fidelity']:.2f} < 0.5(纠缠判据线):该态已被噪声"
                "支配,图中散乱柱体即噪声本身(诚实呈现的负结果"
                + (";读出矩阵近奇异时缓解放大噪声,见缓解失效边界一节" if suffix
                   else "") + ")")
        fig.text(0.5, 0.015, note, ha="center", fontsize=8.5,
                 color=st.ROLE["bad"])
    st.save(fig, f"tomo_{name}_{chip}{suffix}", formats=("png", "svg", "pdf"))
    re_, im_ = np.array(m["rho_re"]), np.array(m["rho_im"])
    with (st.outdir() / f"tomo_{name}_{chip}{suffix}.csv").open("w", newline="") as f:
        w = csv.writer(f)
        n = re_.shape[0]
        labels = [format(i, f"0{int(np.log2(n))}b") for i in range(n)]
        w.writerow(["row\\col", *labels, *(f"im_{x}" for x in labels)])
        for i, lab in enumerate(labels):
            w.writerow([lab, *re_[i].tolist(), *im_[i].tolist()])


def main():
    st.apply_style()
    found = []
    for chip in CHIPS:
        s = latest(chip)
        if not s:
            print(f"跳过 {chip}(无层析数据)")
            continue
        for name, m in s.get("states", {}).items():
            draw_state(chip, name, m)
            if m.get("mitigated"):
                mm = {**m["mitigated"], "qubits": m.get("qubits")}
                draw_state(chip, name, mm, suffix="_mitigated")
        if "bell" in s.get("states", {}):
            found.append((chip, s["states"]["bell"]))
    if len(found) >= 2:                        # 三机对比:每机 Re+Im 一列
        ncol = len(found)
        fig = plt.figure(figsize=(4.8 * ncol, 8.2))
        for i, (chip, m) in enumerate(found):
            axr = fig.add_subplot(2, ncol, i + 1, projection="3d")
            axi = fig.add_subplot(2, ncol, ncol + i + 1, projection="3d")
            plot_state_city(rho_of(m), ax_real=axr, ax_imag=axi)
            axr.set_title(f"{chip} · Re ρ\n{head(m)}", fontsize=9.5,
                          color=st.INK)
            axi.set_title(f"{chip} · Im ρ", fontsize=9.5, color=st.INK)
        fig.subplots_adjust(left=0.0, right=1.0, wspace=0.04, hspace=0.12,
                            top=0.9, bottom=0.03)
        fig.suptitle("同一 Bell 态 · 三台机器断层成像(当日自测)",
                     fontsize=13.5, y=0.97)
        st.save(fig, "tomo_bell_combo", formats=("png", "svg", "pdf"))


if __name__ == "__main__":
    main()
