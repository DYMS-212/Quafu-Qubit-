"""每日数据导出:全部实验的 summary → 整洁 CSV,一天一个文件夹(零额度)。

产出 data_export/<YYYY-MM-DD>/:
  e04_syndrome.csv, e05_repcode.csv, ... e22_variance_panel.csv   (每轮一行/多行,带 run 列)
  e14_chip_eta.csv   每次体检的全芯片逐比特 η/F0/F1(天气图原始数据)
  e14_probe_log.csv  当日体检汇总行
  e02_selfcal.csv    逐比特自校准+vendor对照
  e00_survey.csv / e01_chip.csv / e03_bell.csv / e17_variance.csv
  MANIFEST.md        文件清单+行数(讨论时先看这个)

原始 counts 仍在 results/(唯一真相);CSV 是讨论/复盘用的整洁视图,可随时重生成。

运行: ./qrun tools/export_daily.py            # 全部日期
      ./qrun tools/export_daily.py <日期>     # 只导某天
"""
from __future__ import annotations

import csv
import json
import re
import sys
import time
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))

from qfab import analysis  # noqa: E402

R = PROJ / "results"
OUT = PROJ / "data_export"
DATE_RE = re.compile(r"(\d{4}-\d{2}-\d{2})(?:_(\d{4}))?")
ONLY = sys.argv[1] if len(sys.argv) > 1 else None


def cell(v):
    if isinstance(v, (list, dict)):
        return json.dumps(v, ensure_ascii=False, separators=(",", ":"))
    return v


def write_csv(date: str, name: str, records: list[dict]) -> int:
    if not records:
        return 0
    d = OUT / date
    d.mkdir(parents=True, exist_ok=True)
    cols: list[str] = []
    for r in records:
        for k in r:
            if k not in cols:
                cols.append(k)
    with (d / f"{name}.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in records:
            w.writerow({k: cell(v) for k, v in r.items()})
    return len(records)


def dated_files(exp: str, pattern: str):
    """[(date, run, path)];run=HHMM 或 date 本身。"""
    d = R / exp
    out = []
    if d.exists():
        for f in sorted(d.glob(pattern)):
            m = DATE_RE.search(f.name)
            if not m:
                continue
            date, hm = m.group(1), m.group(2) or ""
            if ONLY and date != ONLY:
                continue
            out.append((date, f"{date}_{hm}" if hm else date, f))
    return out


def rows_of(data):
    return data.get("rows", data) if isinstance(data, dict) else data


# ------------------------------------------------------------- 通用 rows 型
GENERIC = [
    ("e04_syndrome", "summary_*.json"), ("e05_repcode", "summary_*.json"),
    ("e06_postselect", "summary_*.json"), ("e07_flag", "summary_*.json"),
    ("e08_tiling", "summary_*.json"), ("e09_ghz", "summary_*.json"),
    ("e10_code422", "summary_*.json"), ("e11_memory", "summary_*.json"),
    ("e12_decoder", "summary_*.json"), ("e13_mitigation", "summary_*.json"),
    ("e18_t2swap", "summary_*.json"), ("e19_gates", "summary_*.json"),
    ("e20_anchor_ab", "summary_*.json"), ("e21_selection", "summary_*.json"),
    ("e25_iceberg", "summary_*.json"), ("e26_code832", "summary_*.json"),
    ("e22_variance_panel", "summary_*.json"), ("e03_bell", "atlas_Baihua_*.json"),
    ("e24_logical_bell", "summary_*.json"),
]
# 跨机边图谱单列(防混机):文件名=e03_bell_<机>.csv
CROSS_ATLAS = [("e03_bell_Dongling", "e03_bell", "atlas_Dongling_*.json")]


def export_generic(counts: dict):
    per_day: dict[str, dict[str, list]] = {}
    for exp, pattern in GENERIC:
        for date, run, f in dated_files(exp, pattern):
            data = json.loads(f.read_text())
            rows = rows_of(data)
            if not isinstance(rows, list):
                continue
            recs = [{"run": run, **r} for r in rows if isinstance(r, dict)]
            if exp == "e19_gates" and isinstance(data, dict) and data.get("slopes"):
                sl = data["slopes"]
                extra = [{"run": run, "arm": k, **v} for k, v in sl.items()
                         if isinstance(v, dict)]
                per_day.setdefault(date, {}).setdefault("e19_slopes", []).extend(extra)
            per_day.setdefault(date, {}).setdefault(exp, []).extend(recs)
    for date, exps in per_day.items():
        for exp, recs in exps.items():
            counts.setdefault(date, {})[exp] = write_csv(date, exp, recs)


# ------------------------------------------------------------- 特例
def export_cross_atlas(counts: dict):
    for name, exp, pattern in CROSS_ATLAS:
        per_day: dict[str, list] = {}
        for date, run, f in dated_files(exp, pattern):
            rows = rows_of(json.loads(f.read_text()))
            if isinstance(rows, list):
                per_day.setdefault(date, []).extend(
                    {"run": run, **r} for r in rows if isinstance(r, dict))
        for date, recs in per_day.items():
            counts.setdefault(date, {})[name] = write_csv(date, name, recs)


def export_tomography(counts: dict):
    """e23:每 (chip, state) 一行(rho 矩阵不进 CSV,留在 summary json)。"""
    per_day: dict[str, list] = {}
    for date, run, f in dated_files("e23_tomography", "summary_*.json"):
        d = json.loads(f.read_text())
        for state, m in (d.get("states") or {}).items():
            row = {"run": run, "chip": d.get("chip"), "state": state,
                   **{k: m.get(k) for k in ("qubits", "settings", "shots",
                                            "fidelity", "purity", "negativity",
                                            "chsh_max", "mermin")}}
            mit = m.get("mitigated") or {}
            row["fidelity_mitigated"] = mit.get("fidelity")
            per_day.setdefault(date, []).append(row)
    for date, recs in per_day.items():
        counts.setdefault(date, {})["e23_tomography"] = write_csv(
            date, "e23_tomography", recs)


def export_specials(counts: dict):
    # e00 普查:每机一行
    for date, run, f in dated_files("e00_survey", "survey_*.json"):
        d = json.loads(f.read_text())
        recs = [{"run": run, "chip": c, **(info if isinstance(info, dict) else {"status": info})}
                for c, info in d.get("chips", {}).items()]
        counts.setdefault(date, {})["e00_survey"] = write_csv(date, "e00_survey", recs)
    # e01 芯片画像:单行
    for date, run, f in dated_files("e01_chip", "profile_*.json"):
        d = json.loads(f.read_text())
        counts.setdefault(date, {})["e01_chip"] = write_csv(date, "e01_chip", [{"run": run, **d}])
    # e02 自校准:逐比特
    for date, run, f in dated_files("e02_selfcal", "assignment_*.json"):
        recs = [{"run": run, **r} for r in json.loads(f.read_text())]
        counts.setdefault(date, {})["e02_selfcal"] = write_csv(date, "e02_selfcal", recs)
    # e17 方差:bell/rep 各一行(同日多轮合并)
    e17_by_day: dict[str, list] = {}
    for date, run, f in dated_files("e17_variance", "summary_*.json"):
        d = json.loads(f.read_text())
        for kind in ("bell", "rep"):
            v = d.get(kind)
            if isinstance(v, dict):
                e17_by_day.setdefault(date, []).append(
                    {"run": run, "kind": kind, "edge": d.get("edge"),
                     "chain": d.get("chain"), "R": d.get("R"), "shots": d.get("shots"),
                     **{k: v[k] for k in ("mean", "job_sd", "shot_sd_theory",
                                          "excess_ratio") if k in v}})
    for date, recs in e17_by_day.items():
        counts.setdefault(date, {})["e17_variance"] = write_csv(date, "e17_variance", recs)
    # e14:probe_log 行 + 逐比特 η
    log = R / "e14_probe" / "probe_log.jsonl"
    if log.exists():
        by_day: dict[str, list] = {}
        for line in log.read_text().strip().splitlines():
            d = json.loads(line)
            if ONLY and d.get("date") != ONLY:
                continue
            flat = {}
            for k, v in d.items():
                if isinstance(v, dict):
                    for kk, vv in v.items():
                        flat[f"{k}_{kk}"] = vv
                else:
                    flat[k] = v
            by_day.setdefault(d.get("date", "unknown"), []).append(flat)
        for date, recs in by_day.items():
            counts.setdefault(date, {})["e14_probe_log"] = write_csv(
                date, "e14_probe_log", recs)
    probe = R / "e14_probe"
    for day in sorted(x for x in probe.iterdir() if x.is_dir()) if probe.exists() else []:
        if ONLY and day.name != ONLY:
            continue
        recs = []
        f0s = sorted(day.glob("selfcal_prep0__*.json"))
        f1s = sorted(day.glob("selfcal_prep1__*.json"))
        for f0, f1 in zip(f0s, f1s):
            r0, r1 = json.loads(f0.read_text()), json.loads(f1.read_text())
            if not (r0.get("result", {}).get("count") and r1.get("result", {}).get("count")):
                continue
            qubits = r0["meta"]["qubits"]
            asg = analysis.assignment_from_p1(
                analysis.bit_p1(r0["result"]["count"], len(qubits)),
                analysis.bit_p1(r1["result"]["count"], len(qubits)))
            ts = (r0.get("submitted_at") or "")[:16]
            recs += [{"probe_time": ts, "qubit": q, **asg[i]}
                     for i, q in enumerate(qubits)]
        counts.setdefault(day.name, {})["e14_chip_eta"] = write_csv(
            day.name, "e14_chip_eta", recs)


def main():
    counts: dict[str, dict[str, int]] = {}
    export_generic(counts)
    export_cross_atlas(counts)
    export_tomography(counts)
    export_specials(counts)
    stamp = time.strftime("%Y-%m-%d %H:%M")
    for date in sorted(counts):
        lines = [f"# 数据导出 · {date}", "",
                 f"生成于 {stamp};原始 counts 见 results/(唯一真相),本目录可随时重生成。",
                 "", "| 文件 | 行数 |", "|---|---|"]
        for exp in sorted(counts[date]):
            n = counts[date][exp]
            if n:
                lines.append(f"| {exp}.csv | {n} |")
        (OUT / date / "MANIFEST.md").write_text("\n".join(lines) + "\n")
        total = sum(v for v in counts[date].values())
        print(f"{date}: {sum(1 for v in counts[date].values() if v)} 个 CSV,{total} 行")
    print(f"\n已写 {OUT}/")


if __name__ == "__main__":
    main()
