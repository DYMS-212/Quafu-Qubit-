"""e07 · flag-qubit hook 错误检测(高级 FT 原语)。

ZZ 检验中,ancilla 在两次 CX(anc→flag) 之间发生的 X 故障(hook)会同时
污染 syndrome 并传播——flag 比特专职报警:
  时序: CX(D0→anc) · CX(anc→flag) · [此处注错] · CX(anc→flag) · CX(D1→anc)
  预期: X_anc(hook) → flag=1;X_D0/X_D1(普通数据错) → flag=0 且 syndrome 正常报。

选址:只用实测过判别力的站(tools/probe_flag_sites.py 全芯片探针)。
vendor 边保真度与实测判别力不符——旧站(data[20,22]/anc21/flag34)vendor 边
0.978 但实测 flag 本底 0.63、Δ 为负;vendor 分最高的候选站实测反而最差。
故 SITES 表按探针实测排序,主站不可用时顺延,末路回退旧规则并在 selection 标注。

运行: ./qrun experiments/e07_flag.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, chip as chipmod, client, provenance  # noqa: E402
from qfab.qasm import QASM  # noqa: E402

CHIP = "Baihua"
TODAY = time.strftime("%Y-%m-%d")
RUN = time.strftime("%Y-%m-%d_%H%M")   # 同日多轮互不覆盖
EXP = f"e07_flag/{RUN}"
OUT = Path(__file__).resolve().parent.parent / "results" / "e07_flag"
SHOTS = 4096
DATA, ANC = [20, 22], 21          # 旧规则站,仅作末路回退
# 2026-07-21 探针实测可用站(本底<0.15 且 Δ 显著>0),按综合质量排序:
SITES = [dict(data=[72, 58], anc=71, flag=70,   # 本底0.086 Δ+0.813 anc核对0.883
              probe="2026-07-21 Δ=+0.813"),
         dict(data=[85, 97], anc=84, flag=83,   # 本底0.122 Δ+0.707 anc核对0.850
              probe="2026-07-21 Δ=+0.707")]


def circ(case: str) -> str:
    """q0,q1=data; q2=anc; q3=flag。"""
    qc = QASM(4)
    if case == "X_D0":
        qc.x(0)
    if case == "X_D1":
        qc.x(1)
    qc.barrier()
    qc.cx(0, 2)
    qc.cx(2, 3)
    qc.barrier()
    if case == "X_anc":
        qc.x(2)
    qc.barrier()
    qc.cx(2, 3)
    qc.cx(1, 2)
    qc.barrier().measure_all()
    return qc.qasm()


def main():
    print(f"后端状态: {client.status()}")
    ch = chipmod.load(CHIP, refresh=True)
    av = set(ch.available_qubits())
    site = next((s for s in SITES
                 if {*s["data"], s["anc"], s["flag"]} <= av
                 and all((ch.edge_fid(s["anc"], x) or 0) > 0
                         for x in (*s["data"], s["flag"]))), None)
    if site is not None:
        data, anc, flag = site["data"], site["anc"], site["flag"]
        selection = f"probed_site ({site['probe']})"
    else:
        # 末路回退:旧规则(anc21 第三邻居)——该站实测无判别力,数据仅存档
        data, anc = DATA, ANC
        flags = [x for x in ch.neighbors(anc)
                 if x in av and x not in data and (ch.edge_fid(anc, x) or 0) > 0.5]
        if not flags:
            print(f"⚠️ 探针站全不可用且 anc Q{anc} 无第三邻居,中止")
            return
        flag = flags[0]
        selection = "fallback_legacy_site (实测无判别力,谨慎解读)"
    tq = [*data, anc, flag]
    print(f"布局: data={data} anc=Q{anc} flag=Q{flag} [{selection}]")

    cases = ("none", "X_anc", "X_D0", "X_D1")
    specs = [dict(circuit=circ(cs), chip=CHIP, shots=SHOTS,
                  name=f"flag_{cs}", exp=EXP, label=cs,
                  options={"target_qubits": tq},
                  meta={"case": cs, "flag": flag})
             for cs in cases]
    recs = client.run_batch(specs, timeout=1200)

    rows = []
    print("\n===== flag 查找表(期望: X_anc→flag=1;数据错→flag=0)=====")
    for r in recs:
        c = client.counts(r)
        if not c:
            continue
        p1 = analysis.bit_p1(c, 4)
        tot = analysis.shots_of(c)
        f_lo, f_hi = analysis.wilson_ci(round(p1[3] * tot), tot)
        a_lo, a_hi = analysis.wilson_ci(round(p1[2] * tot), tot)
        rows.append({"case": r["meta"]["case"],
                     "p_flag_1": round(p1[3], 4),
                     "p_flag_1_ci": [round(f_lo, 4), round(f_hi, 4)],
                     "p_anc_1": round(p1[2], 4),
                     "p_anc_1_ci": [round(a_lo, 4), round(a_hi, 4)],
                     "shots": tot})
        print(f"  {r['meta']['case']:6s}: P(flag=1)={p1[3]:.4f}  P(anc=1)={p1[2]:.4f}")
    prov = provenance.stamp(
        qubits=tq, tids=provenance.tids_of(recs), chip=CHIP,
        extra={"selection": selection,
               "data": data, "anc": anc, "flag": flag,
               "shots_per_task": SHOTS,
               "ci": "Wilson 95% on P(flag=1) 与 P(anc=1)"})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"summary_{RUN}.json").write_text(
        json.dumps({"provenance": prov, "rows": rows}, ensure_ascii=False, indent=2))
    print(f"已写 {OUT / f'summary_{RUN}.json'}")


if __name__ == "__main__":
    main()
