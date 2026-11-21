"""e13 · 误差缓解全景(赛道名的另一半):
  A. CZ Pauli twirling —— 对症 e03 发现的 ZZ/XX 不对称(相干相位误差→随机化)
  B. ZNE 零噪声外推 —— CZ 折叠 1×/3×/5×,线性外推
  C. vendor correct / DD(XY4/CPMG)—— 平台旋钮对照(在 repcode 锚点上)
  D. 自校准读出纠正 —— 离线,见 e02/e12(此处不重复)

运行: ./qrun experiments/e13_mitigation.py
"""
from __future__ import annotations

import json
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, anchors as anchmod, client, provenance, repcode  # noqa: E402
from qfab.qasm import QASM  # noqa: E402

CHIP = "Baihua"
TODAY = time.strftime("%Y-%m-%d")
RUN = time.strftime("%Y-%m-%d_%H%M")   # 同日多轮互不覆盖
EXP = f"e13_mitigation/{RUN}"
PROJ = Path(__file__).resolve().parent.parent
OUT = PROJ / "results" / "e13_mitigation"
SHOTS = 2048
N_TWIRL = 4
PLAIN_FOLDS = (3, 5, 7)         # plain-ZNE 折叠(+fold=7 检验线性假设)
TWFOLDS = (1, 3, 5, 7)          # twirled-ZNE(物理正确版)折叠
TWFOLD_SEEDS = 4                # 每折叠随机实例数(≥4 才能 bootstrap 出截距 CI)
_A = anchmod.latest(CHIP)        # 选址改用当日锚点(不再硬编码边/链)
BELL_EDGES = [tuple(e) for e in _A.get("bell_edges") or [[67, 68], [73, 74]]]
REP_CHAIN = (_A.get("chains") or {}).get("3") or [8, 7, 6, 5, 4]


def bell_cz(basis: str, twirl_seed: int | None = None, fold: int = 1) -> str:
    """Bell 态(显式 CZ 版,便于 twirl/fold):h0·h1·CZ^fold·h1。"""
    qc = QASM(2)
    qc.h(0).h(1)
    rng = random.Random(twirl_seed) if twirl_seed is not None else None
    for _ in range(fold):
        if rng is not None:
            qc.twirl_cz(0, 1, rng)
        else:
            qc.cz(0, 1)
    qc.h(1)
    if basis == "XX":
        qc.h(0).h(1)
    qc.barrier().measure_all()
    return qc.qasm()


def parity(count: dict) -> float:
    tot = analysis.shots_of(count)
    even = sum(v for k, v in count.items() if k.zfill(2).count("1") % 2 == 0)
    return 2 * even / tot - 1


def zne_intercept_ci(fold_vals: dict, B: int = 1000, seed: int = 0) -> dict:
    """fold_vals={fold:[各seed的parity]}。线性外推到 fold=0 的截距 + bootstrap(重采样seed)CI。
    截距 CI 与 raw 的 CI 是否重叠 → ZNE 有无显著增益(钉死负结果的关键)。"""
    folds = sorted(fold_vals)
    if len(folds) < 2 or any(not fold_vals[f] for f in folds):
        return {}
    rng = random.Random(seed)

    def fit(sample):
        ys = [sum(sample[f]) / len(sample[f]) for f in folds]
        xb, yb = sum(folds) / len(folds), sum(ys) / len(ys)
        denom = sum((x - xb) ** 2 for x in folds) or 1e-9
        k = sum((x - xb) * (y - yb) for x, y in zip(folds, ys)) / denom
        return yb - k * xb

    boots = sorted(fit({f: [rng.choice(fold_vals[f]) for _ in fold_vals[f]] for f in folds})
                   for _ in range(B))
    return {"intercept": round(fit(fold_vals), 4),
            "lo": round(boots[int(B * 0.025)], 4),
            "hi": round(boots[int(B * 0.975)], 4), "n_folds": len(folds)}


def main():
    print(f"后端状态: {client.status()}")
    specs = []
    # A/B: twirling + ZNE(每条 Bell 边 × 两基)
    for (a, b) in BELL_EDGES:
        for basis in ("ZZ", "XX"):
            base = dict(chip=CHIP, shots=SHOTS, exp=EXP,
                        options={"target_qubits": [a, b]})
            specs.append(dict(**base, circuit=bell_cz(basis),
                              name=f"m_raw_{a}_{b}_{basis}",
                              label=f"raw_{a}-{b}_{basis}",
                              meta={"kind": "raw", "edge": [a, b], "basis": basis}))
            for s in range(N_TWIRL):
                specs.append(dict(**base, circuit=bell_cz(basis, twirl_seed=s),
                                  name=f"m_tw{s}_{a}_{b}_{basis}",
                                  label=f"twirl{s}_{a}-{b}_{basis}",
                                  meta={"kind": "twirl", "seed": s,
                                        "edge": [a, b], "basis": basis}))
            for fold in PLAIN_FOLDS:
                specs.append(dict(**base, circuit=bell_cz(basis, fold=fold),
                                  name=f"m_f{fold}_{a}_{b}_{basis}",
                                  label=f"fold{fold}_{a}-{b}_{basis}",
                                  meta={"kind": "fold", "fold": fold,
                                        "edge": [a, b], "basis": basis}))
            # twirled-ZNE(物理正确版):相干误差先被 twirl 随机化,折叠标度才成立
            for fold in TWFOLDS:
                for sd in range(TWFOLD_SEEDS):
                    specs.append(dict(**base,
                                      circuit=bell_cz(basis, twirl_seed=100 * fold + sd,
                                                      fold=fold),
                                      name=f"m_tf{fold}s{sd}_{a}_{b}_{basis}",
                                      label=f"twfold{fold}s{sd}_{a}-{b}_{basis}",
                                      meta={"kind": "twfold", "fold": fold, "seed": sd,
                                            "edge": [a, b], "basis": basis}))
    # C: 平台旋钮(repcode d=3 锚点)
    knobs = {"raw": {}, "vcorrect": {"correct": True},
             "xy4": {"open_dd": "XY4"}, "cpmg": {"open_dd": "CPMG"}}
    for kn, opt in knobs.items():
        for bit in (0, 1):
            specs.append(dict(circuit=repcode.build_t1(3, bit), chip=CHIP,
                              shots=4096, exp=EXP,
                              name=f"m_knob_{kn}_L{bit}",
                              label=f"knob_{kn}_L{bit}",
                              options={"target_qubits": REP_CHAIN, **opt},
                              meta={"kind": "knob", "knob": kn, "bit": bit}))
    print(f"共 {len(specs)} 条线路")
    recs = client.run_batch(specs, timeout=1200)

    by_label = {r["label"]: r for r in recs if client.counts(r)}
    rows = []
    print("\n===== A/B: twirling 与 ZNE(Bell parity)=====")
    for (a, b) in BELL_EDGES:
        for basis in ("ZZ", "XX"):
            raw = by_label.get(f"raw_{a}-{b}_{basis}")
            if not raw:
                continue
            p_raw = parity(client.counts(raw))
            tw = [parity(client.counts(by_label[f"twirl{s}_{a}-{b}_{basis}"]))
                  for s in range(N_TWIRL)
                  if f"twirl{s}_{a}-{b}_{basis}" in by_label]
            p_tw = sum(tw) / len(tw) if tw else None
            folds = {1: p_raw}
            for f in (3, 5):
                lbl = f"fold{f}_{a}-{b}_{basis}"
                if lbl in by_label:
                    folds[f] = parity(client.counts(by_label[lbl]))
            # 线性外推到 fold→0
            zne = None
            if len(folds) == 3:
                xs, ys = list(folds), [folds[x] for x in folds]
                xb, yb = sum(xs) / 3, sum(ys) / 3
                k = sum((x - xb) * (y - yb) for x, y in zip(xs, ys)) / \
                    sum((x - xb) ** 2 for x in xs)
                zne = yb - k * xb
            # twirled-ZNE(物理正确版):每折叠多 seed → 截距 + bootstrap CI
            tw_fold_vals = {}
            for fold in TWFOLDS:
                vals_f = [parity(client.counts(by_label[f"twfold{fold}s{sd}_{a}-{b}_{basis}"]))
                          for sd in range(TWFOLD_SEEDS)
                          if f"twfold{fold}s{sd}_{a}-{b}_{basis}" in by_label]
                if vals_f:
                    tw_fold_vals[fold] = vals_f
            tw_zne = zne_intercept_ci(tw_fold_vals)
            raw_c = client.counts(raw)
            even = sum(v for k, v in raw_c.items() if k.zfill(2).count("1") % 2 == 0)
            lo, hi = analysis.wilson_ci(even, analysis.shots_of(raw_c))
            raw_ci_par = [round(2 * lo - 1, 4), round(2 * hi - 1, 4)]  # 映射到 parity 尺度
            overlap = bool(tw_zne) and tw_zne["lo"] <= raw_ci_par[1] and tw_zne["hi"] >= raw_ci_par[0]
            rows.append({"edge": [a, b], "basis": basis,
                         "raw": round(p_raw, 4), "raw_ci": raw_ci_par,
                         "twirl": round(p_tw, 4) if p_tw is not None else None,
                         "folds_plain": {str(k): round(v, 4) for k, v in folds.items()},
                         "zne_plain": round(zne, 4) if zne is not None else None,
                         "tw_folds_mean": {str(k): round(sum(v) / len(v), 4)
                                           for k, v in tw_fold_vals.items()},
                         "tw_zne": tw_zne,
                         "zne_significant": (not overlap) if tw_zne else None})
            print(f"  {a}-{b} {basis}: raw={p_raw:.4f}{raw_ci_par} | "
                  f"tw_zne截距={tw_zne.get('intercept') if tw_zne else None}"
                  f"[{tw_zne.get('lo') if tw_zne else '?'},{tw_zne.get('hi') if tw_zne else '?'}]"
                  f" → 与raw {'重叠(无显著增益)' if overlap else '不重叠(有增益)'}")
    print("\n===== C: 平台旋钮(repcode d=3 多数表决错误率)=====")
    for kn in ("raw", "vcorrect", "xy4", "cpmg"):
        vals = {}
        for bit in (0, 1):
            r = by_label.get(f"knob_{kn}_L{bit}")
            if r:
                c = client.counts(r)
                # vendor correct 时平台可能回 corrected 字段,统一用 count
                vals[f"L{bit}"] = round(
                    repcode.logical_error_majority(c, 3, bit), 5)
        rows.append({"kind": "knob", "knob": kn, **vals})
        print(f"  {kn:9s}: {vals}")
    OUT.mkdir(parents=True, exist_ok=True)
    summary = {"provenance": provenance.stamp(
                   anchor=_A, chip=CHIP, tids=provenance.tids_of(recs),
                   qubits={"bell_edges": [list(e) for e in BELL_EDGES], "rep_chain": REP_CHAIN},
                   extra={"twfold_seeds": TWFOLD_SEEDS, "twfolds": list(TWFOLDS),
                          "note": "twirled-ZNE=物理正确版(截距带bootstrap CI);plain-ZNE 作相干误差下失效对照"}),
               "rows": rows}
    (OUT / f"summary_{RUN}.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"已写 {OUT / f'summary_{RUN}.json'}")


if __name__ == "__main__":
    main()
