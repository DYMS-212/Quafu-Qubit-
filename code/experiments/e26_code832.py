"""e26 · [[8,3,2]] 立方码:3 逻辑比特寄存器 + transversal T = 逻辑 CCZ(magic 门探针)。

两条臂(全部 statevector 现场自校准,绝不手推):

**寄存器臂**:|0̄0̄0̄⟩ = GHZ₈(路径中点向两侧,深度 4),8 个逻辑基态由
X̄_b(立方体面上的物理 X 层)写入;全 Z 测量 → 4 重独立稳定子奇偶后选
(Z-面 ×3 + Z⊗8)→ 读出 3 个逻辑比特。8 物理比特存 3 逻辑比特、
四重检错——比 [[4,2,2]] 高一档的检错寄存器。

**CCZ 臂(magic)**:|+̄+̄+̄⟩(RM(1,3) 均匀叠加,CSS 主元法制备 + 路径 SWAP
布线)→ 按顶点奇偶打 T/T†(单比特层!)= 逻辑 CCZ → 三个混合基设置
(面 v_b=1 测 X、其余测 Z)读 hypergraph 稳定子 K_b = X̄_b·CZ̄:
相关子 E[x̄_b·(-1)^{z̄'z̄''}] 理想 = ±1;同批无 T 对照臂给出零假设参考。
CCZ 非 Clifford——"一层单比特门实现三逻辑比特 magic 门"是本码的招牌;
若门质量不够、相关子被制备深度吃掉,如实报边界(制备深度见 provenance)。

立方体顶点 v∈{0..7}(3 bit),Gray 序放到路径上(面尽量连续,SWAP 最少)。
运行: ./qrun experiments/e26_code832.py [chip]      # 默认 Baihua
      ./qrun experiments/e26_code832.py selftest    # 零额度 statevector 自检
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, anchors as anchmod, client, provenance  # noqa: E402
from qfab.qasm import QASM  # noqa: E402

CHIP = sys.argv[1] if len(sys.argv) > 1 else "Baihua"
RUN = time.strftime("%Y-%m-%d_%H%M")
TAG = RUN if CHIP == "Baihua" else f"{CHIP}_{RUN}"
EXP = f"e26_code832/{TAG}"
OUT = Path(__file__).resolve().parent.parent / "results" / "e26_code832"
SHOTS = 4096
PI = 3.141592653589793

GRAY = [0, 1, 3, 2, 6, 7, 5, 4]                  # 路径位 i ← 顶点 GRAY[i]
POS = {v: i for i, v in enumerate(GRAY)}          # 顶点 → 路径位


def face(b: int) -> list[int]:
    """面 {v: v_b=1} 的路径位集合(= X̄_b 支撑;补面 = Z 检查之一)。"""
    return sorted(POS[v] for v in range(8) if (v >> b) & 1)


Z_CHECKS = [face(0), face(1), face(2), list(range(8))]   # 4 独立 Z 稳定子
ZBAR = [(POS[0], POS[1 << b]) for b in range(3)]          # Z̄_b = Z_000 Z_e_b


# ------------------------------------------------------------------ 线路
def ghz8(qc: QASM):
    qc.h(4)
    for i in range(4, 7):
        qc.cx(i, i + 1)
    for i in range(4, 0, -1):
        qc.cx(i, i - 1)


def xbar_pattern(x: int, y: int, z: int) -> set[int]:
    s: set[int] = set()
    for b, on in enumerate((x, y, z)):
        if on:
            s ^= set(face(b))
    return s


def circ_register(x: int, y: int, z: int) -> str:
    qc = QASM(8)
    ghz8(qc)
    for q in sorted(xbar_pattern(x, y, z)):
        qc.x(q)
    qc.barrier().measure_all()
    return qc.qasm()


def _row_spread(qc: QASM, p: int, supp: list[int]) -> int:
    """主元 p 的行展开:主元内容沿路径单程右走(SWAP),路过每个目标就近 CX,
    最后原路走回(恢复被临时挪位的比特)。返回发射的 CX 数。"""
    def swap(x, y):
        qc.cx(x, y).cx(y, x).cx(x, y)
    cur, n = p, 0
    for j in sorted(supp):
        while cur < j - 1:
            swap(cur, cur + 1)
            cur += 1
            n += 3
        qc.cx(cur, j)
        n += 1
    while cur > p:
        swap(cur - 1, cur)
        cur -= 1
        n += 3
    return n


def circ_plus(with_t: bool, t_sign: int, setting: int) -> tuple[str, int]:
    """|+̄+̄+̄⟩(CSS 主元法 RREF + 路径布线)→(可选)T 层 → setting b 混合基测量。
    返回 (qasm, 双比特门数)。t_sign=+1: 偶奇偶顶点打 T;-1: 反过来。
    setting=-1 表示全 Z 测量(自检用)。"""
    # RREF(RM(1,3),Gray 路径序):主元 + 其余支撑
    # (selftest 以 Z 基支撑 16 词均匀 + CCZ 相关子确定性核验)
    rows = [(0, [3, 5, 6]), (1, [3, 5, 7]), (2, [3, 6, 7]), (4, [5, 6, 7])]
    qc = QASM(8)
    n2q = 0
    for p, _ in rows:
        qc.h(p)
    for p, supp in rows:
        n2q += _row_spread(qc, p, supp)
    qc.barrier()
    if with_t:
        for v in range(8):
            odd = bin(v).count("1") % 2
            theta = (PI / 4) * (t_sign if odd == 0 else -t_sign)
            qc.rz(theta, POS[v])
        qc.barrier()
    if setting >= 0:
        for q in face(setting):
            qc.h(q)
    qc.barrier().measure_all()
    return qc.qasm(), n2q


# ------------------------------------------------------------------ 分析
def eval_register(count: dict, exp_checks: list[int],
                  exp_log: tuple[int, int, int]) -> dict:
    tot = acc = any_err = 0
    per = [0, 0, 0]
    for key, v in count.items():
        b = analysis.bits_of(key, 8)
        tot += v
        if any(sum(b[q] for q in chk) % 2 != e
               for chk, e in zip(Z_CHECKS, exp_checks)):
            continue
        acc += v
        bad = False
        for i, (qa, qb) in enumerate(ZBAR):
            if (b[qa] ^ b[qb]) != exp_log[i]:
                per[i] += v
                bad = True
        if bad:
            any_err += v
    return {"acceptance": round(acc / (tot or 1), 4), "n_acc": acc,
            "err_any": round(any_err / acc, 6) if acc else None,
            "per_logical": [round(x / acc, 6) if acc else None for x in per]}


def eval_ccz(count: dict, setting: int, exp_check: int) -> dict:
    """setting b:面 F_b 已转 X 基。相关子 = x̄_b · (-1)^{z̄' z̄''};
    后选:补面(v_b=0)的 Z 奇偶 == 理想。"""
    fb = set(face(setting))
    comp = [q for q in range(8) if q not in fb]
    others = [b_ for b_ in range(3) if b_ != setting]
    tot = acc = 0
    s = 0
    for key, v in count.items():
        b = analysis.bits_of(key, 8)
        tot += v
        if sum(b[q] for q in comp) % 2 != exp_check:
            continue
        acc += v
        xbar = sum(b[q] for q in fb) % 2
        z1 = b[ZBAR[others[0]][0]] ^ b[ZBAR[others[0]][1]]
        z2 = b[ZBAR[others[1]][0]] ^ b[ZBAR[others[1]][1]]
        # K_b = X̄_b·CZ̄:相关子 = (-1)^{x̄_b} · (-1)^{z̄'∧z̄''}
        s += v * (1 - 2 * xbar) * (1 - 2 * (z1 & z2))
    return {"acceptance": round(acc / (tot or 1), 4), "n_acc": acc,
            "corr": round(s / acc, 4) if acc else None}


def ideal_register(x, y, z):
    from tools.verify_circuits_offline import simulate
    cnt = simulate(circ_register(x, y, z))
    checks, logs = set(), set()
    for key in cnt:
        b = analysis.bits_of(key, 8)
        checks.add(tuple(sum(b[q] for q in chk) % 2 for chk in Z_CHECKS))
        logs.add(tuple(b[qa] ^ b[qb] for qa, qb in ZBAR))
    assert len(checks) == 1 and len(logs) == 1, (x, y, z, checks, logs)
    return list(checks.pop()), logs.pop()


def ideal_ccz():
    """选定 T 符号模式:要求 3 个设置的相关子全部理想 ±1(确定值);
    返回 (t_sign, {setting: (check_parity, ideal_corr)}, 无T对照理想值)。"""
    from tools.verify_circuits_offline import simulate
    for t_sign in (1, -1):
        cal, ok = {}, True
        for b in range(3):
            qasm, _ = circ_plus(True, t_sign, b)
            cnt = simulate(qasm)
            comp = [q for q in range(8) if q not in set(face(b))]
            pars = {tuple([sum(analysis.bits_of(k, 8)[q] for q in comp) % 2])
                    for k in cnt}
            if len(pars) != 1:
                ok = False
                break
            par = pars.pop()[0]
            m = eval_ccz(cnt, b, par)
            if m["corr"] is None or abs(abs(m["corr"]) - 1) > 1e-9 \
                    or abs(m["acceptance"] - 1) > 1e-9:
                ok = False
                break
            cal[b] = (par, round(m["corr"]))
        if ok:
            ctrl = {}
            for b in range(3):
                qasm, _ = circ_plus(False, t_sign, b)
                cnt = simulate(qasm)
                ctrl[b] = eval_ccz(cnt, b, cal[b][0])["corr"]
            return t_sign, cal, ctrl
    raise AssertionError("两种 T 符号模式都未产生确定的 CCZ 相关子——设计有误")


def selftest():
    for x in (0, 1):
        for y in (0, 1):
            for z in (0, 1):
                chk, lg = ideal_register(x, y, z)
                ok = lg == (x, y, z)
                print(f"  {'✅' if ok else '❌'} 寄存器 |{x}{y}{z}⟩_L: "
                      f"检查奇偶={chk} 读出={lg}")
                if not ok:
                    sys.exit(1)
    from tools.verify_circuits_offline import simulate
    qasm, n2q = circ_plus(False, 1, -1)
    cnt = simulate(qasm)
    ok = len(cnt) == 16 and max(cnt.values()) - min(cnt.values()) < 1e-9
    print(f"  {'✅' if ok else '❌'} |+̄+̄+̄⟩ Z 基支撑 = RM(1,3) 均匀 16 词"
          f"(制备双比特门 {n2q} 个)")
    if not ok:
        sys.exit(1)
    t_sign, cal, ctrl = ideal_ccz()
    print(f"  ✅ CCZ:T 符号模式={'偶+奇-' if t_sign == 1 else '偶-奇+'};"
          f"理想相关子 {dict((b, c[1]) for b, c in cal.items())};"
          f"无T对照 {ctrl}")
    print("selftest 全部通过 ✅")


# ------------------------------------------------------------------ 主流程
def main():
    if CHIP == "selftest":
        selftest()
        return
    st = client.status()
    print(f"后端状态: {st}")
    if not isinstance(st.get(CHIP), int):
        print(f"⚠️ {CHIP} 当前不可用({st.get(CHIP)}),中止。")
        return
    A = anchmod.latest(CHIP)
    chain = max((A.get("chains") or {}).values(), key=len, default=[])
    if len(chain) < 8:
        print("⚠️ 锚点链不足 8 比特,中止。")
        return
    tq = chain[:8]
    t_sign, cal, ctrl_ideal = ideal_ccz()
    reg_cal = {(x, y, z): ideal_register(x, y, z)
               for x in (0, 1) for y in (0, 1) for z in (0, 1)}

    specs = []
    for (x, y, z) in reg_cal:
        specs.append(dict(circuit=circ_register(x, y, z), chip=CHIP, shots=SHOTS,
                          name=f"reg_{x}{y}{z}", exp=EXP, label=f"reg_{x}{y}{z}",
                          options={"target_qubits": tq},
                          meta={"kind": "reg", "state": [x, y, z]}))
    n2q_prep = None
    for with_t in (True, False):
        for b in range(3):
            qasm, n2q = circ_plus(with_t, t_sign, b)
            n2q_prep = n2q
            lab = f"ccz_{'T' if with_t else 'noT'}_s{b}"
            specs.append(dict(circuit=qasm, chip=CHIP, shots=SHOTS,
                              name=lab, exp=EXP, label=lab,
                              options={"target_qubits": tq},
                              meta={"kind": "ccz", "with_t": with_t,
                                    "setting": b}))
    print(f"共 {len(specs)} 条线路(寄存器 8 + CCZ 6;CCZ 制备 2q 门 {n2q_prep})")
    recs = client.run_batch(specs, timeout=1200)

    rows = []
    for r in recs:
        c = client.counts(r)
        if not c:
            continue
        m = r["meta"]
        if m["kind"] == "reg":
            chk, lg = reg_cal[tuple(m["state"])]
            rows.append({**m, "chip": CHIP, "qubits": tq, "tid": r.get("tid"),
                         **eval_register(c, chk, lg)})
        else:
            b = m["setting"]
            out = eval_ccz(c, b, cal[b][0])
            if out["n_acc"]:
                lo, hi = analysis.wilson_ci(
                    round((out["corr"] + 1) / 2 * out["n_acc"]), out["n_acc"])
                out["corr_ci"] = [round(2 * lo - 1, 4), round(2 * hi - 1, 4)]
            rows.append({**m, "chip": CHIP, "qubits": tq, "tid": r.get("tid"),
                         "ideal": cal[b][1] if m["with_t"] else ctrl_ideal[b],
                         **out})

    OUT.mkdir(parents=True, exist_ok=True)
    summary = {"provenance": provenance.stamp(
                   anchor=A, chip=CHIP, tids=provenance.tids_of(recs),
                   qubits={"path": tq},
                   extra={"gray_order": GRAY, "t_sign": t_sign,
                          "ccz_prep_2q_gates": n2q_prep,
                          "note": "寄存器臂后选=4 独立 Z 稳定子;CCZ 臂后选="
                                  "补面 Z 奇偶;全部理想值 statevector 自校准"}),
               "rows": rows}
    (OUT / f"summary_{TAG}.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2))

    reg = [r for r in rows if r["kind"] == "reg"]
    if reg:
        accs = [r["acceptance"] for r in reg]
        errs = [r["err_any"] for r in reg if r["err_any"] is not None]
        print(f"\n寄存器臂(8 态):接受率 {min(accs):.3f}-{max(accs):.3f},"
              f"后选后任一逻辑错 {min(errs):.4f}-{max(errs):.4f}")
    for r in rows:
        if r["kind"] == "ccz":
            print(f"  CCZ {'T ' if r['with_t'] else '对照'} 设置{r['setting']}: "
                  f"corr={r['corr']}(理想 {r['ideal']})acc={r['acceptance']}")
    print(f"已写 {OUT / f'summary_{TAG}.json'}")


if __name__ == "__main__":
    main()
