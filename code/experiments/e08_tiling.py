"""e08 · 并行 tiling + 串扰:K 个 d=3 码同芯片并行 vs 逐个单跑。

选址:e02 最新 η 打分,贪心找 K 条不相交 5 比特链(链内边 vendor CZ>0.5)。
对照:packed(一条线路全并行)vs isolated(每 tile 单独线路),
     各测 L0 无注错 + L0 中位注错,packed−isolated 差 = 串扰。
**全芯片档**:放宽 η 门槛(>0.3)把能铺的 tile 全铺满 → "整片芯片同时跑
纠错"的规模上限 + 逐 tile 错误率地图(质量档 A/B 对照仍用 η>0.5 严选)。

运行: ./qrun experiments/e08_tiling.py [K]   # 默认尽量多
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import analysis, chip as chipmod, client, provenance, repcode  # noqa: E402
from qfab.qasm import QASM  # noqa: E402

CHIP = "Baihua"
TODAY = time.strftime("%Y-%m-%d")
RUN = time.strftime("%Y-%m-%d_%H%M")   # 同日多轮互不覆盖
EXP = f"e08_tiling/{RUN}"
PROJ = Path(__file__).resolve().parent.parent
OUT = PROJ / "results" / "e08_tiling"
SHOTS = 4096
K_MAX = int(sys.argv[1]) if len(sys.argv) > 1 else 99


def load_eta() -> tuple[dict[int, float], str]:
    """返回 (η 表, 来源文件名);来源文件名进 provenance,记清选址依据。"""
    files = sorted((PROJ / "results" / "e02_selfcal").glob("assignment_*.json"))
    rows = json.loads(files[-1].read_text())
    return {r["qubit"]: r["eta"] for r in rows}, files[-1].name


def find_tiles(ch, eta, k_max: int, eta_min: float = 0.5) -> list[list[int]]:
    """贪心:按链 min-η 分数找不相交 5 比特路径。"""
    av = set(ch.available_qubits()) & {q for q, e in eta.items() if e > eta_min}
    adj = {q: [x for x in ch.neighbors(q) if x in av
               and (ch.edge_fid(q, x) or 0) > 0.5] for q in av}
    paths = []
    for start in av:
        stack = [[start]]
        while stack:
            p = stack.pop()
            if len(p) == 5:
                paths.append((min(eta[q] for q in p), p))
                continue
            for x in adj[p[-1]]:
                if x not in p:
                    stack.append(p + [x])
    paths.sort(key=lambda t: -t[0])
    tiles, used = [], set()
    for _, p in paths:
        if not (set(p) & used):
            tiles.append(p)
            used |= set(p)
            if len(tiles) >= k_max:
                break
    return tiles


def packed_circuit(k: int, inject_mid: bool) -> str:
    """k 个 d=3 tile 串在一条线路(每 tile 5 逻辑位)。"""
    qc = QASM(5 * k)
    for t in range(k):
        o = 5 * t
        if inject_mid:
            qc.x(o + 2)
        qc.barrier()
        qc.cx(o + 0, o + 1).cx(o + 2, o + 1)
        qc.cx(o + 2, o + 3).cx(o + 4, o + 3)
    qc.barrier().measure_all()
    return qc.qasm()


def tile_err(count: dict, k: int, tile_idx: int) -> float:
    """从 packed count 抽出第 tile_idx 个 tile 的 data 位做多数表决(L0)。"""
    n = 5 * k
    tot = err = 0
    for key, v in count.items():
        bits = analysis.bits_of(key, n)
        o = 5 * tile_idx
        data = [bits[o], bits[o + 2], bits[o + 4]]
        tot += v
        if repcode.majority(data) != 0:
            err += v
    return err / (tot or 1)


def main():
    print(f"后端状态: {client.status()}")
    ch = chipmod.load(CHIP, refresh=True)
    eta, eta_src = load_eta()
    tiles = find_tiles(ch, eta, K_MAX)
    k = len(tiles)
    print(f"选出 {k} 个不相交 d=3 tile(共 {5 * k} 比特):")
    for t in tiles:
        print(f"  {t} (min η={min(eta[q] for q in t):.3f})")

    # 全芯片档:η>0.3 全铺满(只跑 packed,规模上限 + 逐 tile 地图)
    tiles_full = find_tiles(ch, eta, K_MAX, eta_min=0.3)
    kf = len(tiles_full)
    print(f"全芯片档:{kf} 个 tile(共 {5 * kf} 比特,"
          f"占可用比特 {5 * kf}/{len(ch.available_qubits())})")

    specs = []
    for inj in (False, True):
        tag = "injmid" if inj else "none"
        specs.append(dict(circuit=packed_circuit(k, inj), chip=CHIP, shots=SHOTS,
                          name=f"tile_packed_{tag}", exp=EXP,
                          label=f"packed_{tag}",
                          options={"target_qubits": [q for t in tiles for q in t]},
                          meta={"kind": "packed", "inject": inj, "tiles": tiles}))
    if kf > k:
        specs.append(dict(circuit=packed_circuit(kf, False), chip=CHIP, shots=SHOTS,
                          name="tile_fullchip", exp=EXP, label="fullchip_none",
                          options={"target_qubits": [q for t in tiles_full for q in t]},
                          meta={"kind": "fullchip", "inject": False,
                                "tiles": tiles_full}))
        for ti, t in enumerate(tiles):
            specs.append(dict(circuit=repcode.build_t1(3, 0, [1] if inj else []),
                              chip=CHIP, shots=SHOTS,
                              name=f"tile_iso{ti}_{tag}", exp=EXP,
                              label=f"iso{ti}_{tag}",
                              options={"target_qubits": t},
                              meta={"kind": "iso", "inject": inj, "tile_idx": ti,
                                    "tile": t}))
    print(f"共 {len(specs)} 条线路")
    recs = client.run_batch(specs, timeout=1200)

    by = {r["label"]: r for r in recs if client.counts(r)}
    rows = []
    print("\n===== packed vs isolated(L0 多数表决错误率)=====")
    for inj in (False, True):
        tag = "injmid" if inj else "none"
        pk = by.get(f"packed_{tag}")
        for ti, t in enumerate(tiles):
            iso = by.get(f"iso{ti}_{tag}")
            if not (pk and iso):
                continue
            c_pk, c_iso = client.counts(pk), client.counts(iso)
            n_pk, n_iso = analysis.shots_of(c_pk), analysis.shots_of(c_iso)
            e_pk = tile_err(c_pk, k, ti)
            e_iso = repcode.logical_error_majority(c_iso, 3, 0)
            pk_lo, pk_hi = analysis.wilson_ci(round(e_pk * n_pk), n_pk)
            iso_lo, iso_hi = analysis.wilson_ci(round(e_iso * n_iso), n_iso)
            # 串扰 Δ = 两独立比例之差 → Wald 双比例 95% 区间(工具箱只有单比例 Wilson)
            se = (e_pk * (1 - e_pk) / n_pk + e_iso * (1 - e_iso) / n_iso) ** 0.5
            dv = e_pk - e_iso
            rows.append({"tile": t, "inject": inj,
                         "packed": round(e_pk, 5),
                         "packed_ci": [round(pk_lo, 5), round(pk_hi, 5)],
                         "isolated": round(e_iso, 5),
                         "isolated_ci": [round(iso_lo, 5), round(iso_hi, 5)],
                         "delta": round(dv, 5),
                         "delta_ci": [round(dv - 1.96 * se, 5),
                                      round(dv + 1.96 * se, 5)]})
            print(f"  tile{ti} {tag:7s}: packed={e_pk:.4f} iso={e_iso:.4f} "
                  f"Δ={dv:+.4f} [{dv - 1.96 * se:+.4f}, {dv + 1.96 * se:+.4f}]")
    fc = by.get("fullchip_none")
    if fc:
        c_fc = client.counts(fc)
        errs = [tile_err(c_fc, kf, ti) for ti in range(kf)]
        med = sorted(errs)[kf // 2]
        rows.append({"kind": "fullchip", "n_tiles": kf, "n_qubits": 5 * kf,
                     "tiles": tiles_full,
                     "per_tile_err": [round(e, 5) for e in errs],
                     "median_err": round(med, 5)})
        print(f"  全芯片档: {kf} 个码 / {5 * kf} 比特同时运行,L0 中位 {med:.4f},"
              f"最好 {min(errs):.4f},最差 {max(errs):.4f}")

    # 溯源戳:选址来源(e02 自测 η 的具体文件)+ 贪心 tiling + tile 比特 + 任务 id。
    prov = provenance.stamp(
        anchor=None, qubits=[q for t in tiles for q in t],
        tids=provenance.tids_of(recs), chip=CHIP,
        extra={"site_selection": "e02_selfcal_eta + greedy_disjoint_5q_tiling",
               "eta_source": eta_src, "chip_calibration_time": ch.calibration_time,
               "n_tiles": k, "tiles": tiles,
               "n_tiles_fullchip": kf, "tiles_fullchip": tiles_full})
    out = {"provenance": prov, "rows": rows}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"summary_{RUN}.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2))
    print(f"已写 {OUT / f'summary_{RUN}.json'}")


if __name__ == "__main__":
    main()
