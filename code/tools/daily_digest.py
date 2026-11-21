"""每日摘要:自动生成《半决赛报告.md》附录 A 登记表的当日行(直接复制粘贴)。

运行: ./qrun tools/daily_digest.py [YYYY-MM-DD]   # 省略=今天;可指定日期补历史日
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

PROJ = Path(__file__).resolve().parent.parent
R = PROJ / "results"
# 可选 argv[1]=YYYY-MM-DD 覆盖(补历史日);否则今天。下游所有按日过滤自动跟随。
TODAY = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1][:2] == "20" else time.strftime("%Y-%m-%d")
NA = "🔲"


def latest_summary(exp: str, pattern: str = "summary_2*.json"):
    d = R / exp
    if not d.exists():
        return None
    files = sorted(d.glob(pattern))
    # 只认今日产出:某实验今天崩了 → 返回 None → 该列诚实显示 🔲,绝不用昨天的数冒充今天
    files = [f for f in files if TODAY in f.name]
    if not files:
        return None
    data = json.loads(files[-1].read_text())
    return data.get("rows", data) if isinstance(data, dict) else data  # 兼容 {"provenance","rows"} 新形


def all_today(exp: str, pattern: str = "summary_2*.json"):
    """今天所有轮次的 rows(每个文件一份)——供多轮 best/mean±CI 聚合。"""
    d = R / exp
    out = []
    if d.exists():
        for f in sorted(d.glob(pattern)):
            if TODAY in f.name:
                data = json.loads(f.read_text())
                out.append(data.get("rows", data) if isinstance(data, dict) else data)
    return out


def get(fn, default=NA):
    try:
        v = fn()
        return v if v is not None else default
    except Exception:  # noqa: BLE001
        return default


def main():
    # η中位 + 冻结链 L1(体检)
    def probe():
        lines = (R / "e14_probe" / "probe_log.jsonl").read_text().strip().splitlines()
        days = [json.loads(x) for x in lines if json.loads(x).get("date") == TODAY]
        return days[-1] if days else None
    p = get(probe, None)
    eta = get(lambda: p["selfcal"]["eta_median"]) if p else NA

    # 最优链 d7 L1(e05)
    def d7l1():
        rows = latest_summary("e05_repcode")
        v = [r for r in rows if r["code"] == "bit" and r["d"] == 7 and r["bit"] == 1]
        return v[0]["err"] if v else None
    d7 = get(d7l1)

    # e10 最好星 + 最严基线(fusion 行有 ps_err_perq 但带 kind,须排除)
    def e10_pair():
        rows = latest_summary("e10_code422")
        stars: dict = {}
        for r in rows:
            if "ps_err_perq" in r and "kind" not in r:   # 只 code422 行
                stars.setdefault(r["star"], []).append(r["ps_err_perq"])
        best_star = min((sum(v) / len(v) for v in stars.values()), default=None)
        per_q: dict = {}
        for r in rows:
            if r.get("kind") == "phys":
                per_q.setdefault(r["qubit"], []).append(r["err"])
        best_phys = min((sum(v) / len(v) for v in per_q.values()), default=None)
        return best_star, best_phys
    s10, b10 = get(lambda: e10_pair()[0]), get(lambda: e10_pair()[1])

    # e10 融合制备(FT 验证)最好路径
    def e10_fusion():
        rows = latest_summary("e10_code422")
        paths: dict = {}
        for r in rows:
            if r.get("kind") == "fusion":
                paths.setdefault(r["vq"], []).append(r["ps_err_perq"])
        return min((sum(v) / len(v) for v in paths.values()), default=None)
    f10 = get(e10_fusion)

    # e11 d7 K32 / phys(bit 码;旧数据无 code 字段=bit)
    def e11_pair(code="bit"):
        rows = latest_summary("e11_memory")
        enc = [r for r in rows if r.get("kind") == "enc" and r["d"] == 7
               and r["bit"] == 1 and r["K"] == 32
               and r.get("code", "bit") == code]
        phys = [r["err"] for r in rows if r.get("kind") == "phys"
                and r["bit"] == 1 and r["K"] == 32
                and r.get("code", "bit") == code]
        return (enc[0]["err"] if enc else None), (min(phys) if phys else None)
    e11e, e11p = get(lambda: e11_pair()[0]), get(lambda: e11_pair()[1])
    e11pe, e11pp = get(lambda: e11_pair("phase")[0]), get(lambda: e11_pair("phase")[1])

    # e19 交叉 k(第一个 enc<bare 的 k)+ 斜率比(主指标)
    def cross():
        rows = latest_summary("e19_gates")
        for r in sorted(rows, key=lambda x: x["k"]):
            if r.get("ps_err_perq") is not None and r.get("bare_err_perq") is not None \
                    and r["ps_err_perq"] < r["bare_err_perq"]:
                return r["k"]
        return "未交叉"
    k19 = get(cross)

    def slope_ratio():
        d = R / "e19_gates"
        files = [f for f in sorted(d.glob("summary_*.json")) if TODAY in f.name]
        data = json.loads(files[-1].read_text())
        v = data.get("slopes", {}).get("ratio_bare_over_enc")
        return f"{v:.1f}" if v else None
    r19 = get(slope_ratio)
    k19 = f"{k19}|r×{r19}" if r19 != NA else k19

    # e18 mid 特征(最小 p_expected among *_mid)
    def e18v():
        rows = latest_summary("e18_t2swap")  # latest_summary 已自动拆出 rows 列表
        vals = [r["p_expected"] for r in rows if "_mid" in r["case"]]
        return f"{min(vals):.2f}-{max(vals):.2f}" if vals else None
    v18 = get(e18v)

    # Shenglian rep d3 L1(e15 定点质量套件,当日行)
    def shl():
        f = R / "e15_crossbackend" / "crossbackend_log.jsonl"
        rows_ = [json.loads(x) for x in f.read_text().splitlines() if x.strip()]
        rows_ = [d for d in rows_
                 if d.get("chip") == "Shenglian" and d.get("date") == TODAY]
        return rows_[-1].get("rep_d3_L1") if rows_ else None
    sl = get(shl)

    from qfab import client
    used = get(lambda: client.quota_used_all(TODAY))  # A+B+... 合计(failover 后不漏 B),按指定日

    def fmt(x):
        return f"{x:.4f}" if isinstance(x, float) else str(x)
    row = (f"| {TODAY[5:]} | {fmt(eta)} | {fmt(d7)} | {fmt(s10)} | {fmt(b10)} "
           f"| {fmt(e11e)} / {fmt(e11p)} | {k19} | {v18} | {sl} | 额度 {used} |")
    print("附录 A 当日行(复制粘贴到 半决赛报告.md):\n")
    print(row)
    # e23/e24:三机层析与逻辑 Bell(当日)
    def tomo_bell(chip):
        files = [f for f in sorted((R / "e23_tomography").glob(f"summary_{chip}_*.json"))
                 if TODAY in f.name]
        s = json.loads(files[-1].read_text())["states"].get("bell", {})
        return s.get("fidelity")

    def lbell(chip):
        files = [f for f in sorted((R / "e24_logical_bell").glob(f"summary_{chip}_*.json"))
                 if TODAY in f.name]
        a = json.loads(files[-1].read_text())["arms"]
        return a.get("logical", {}).get("fidelity"), a.get("physical", {}).get("fidelity")

    print(f"\n新臂(备注列素材): e10融合最好路径={fmt(f10)} "
          f"| e11相位 d7K32 enc/phys={fmt(e11pe)}/{fmt(e11pp)}")
    for chip in ("Baihua", "Dongling", "Shenglian"):
        tb = get(lambda c=chip: tomo_bell(c))
        lb = get(lambda c=chip: lbell(c), (NA, NA))
        print(f"  {chip}: 层析BellF={fmt(tb)} | 逻辑Bell F={fmt(lb[0])} vs 物理 {fmt(lb[1])}")

    # 多轮聚合(当日全部轮次)——让"末轮不代表全天"显式暴露;冻结后按此出终版
    from qfab import analysis

    def pool(exp, ex):
        vals = [ex(rows) for rows in all_today(exp)]
        return analysis.pool_ci([v for v in vals if isinstance(v, (int, float))])

    aggs = [
        ("e05 最优链 d7L1", "e05_repcode",
         lambda rows: next((r["err"] for r in rows if r.get("code") == "bit"
                            and r.get("d") == 7 and r.get("bit") == 1), None)),
        ("e11 编码 d7 K32", "e11_memory",
         lambda rows: next((r["err"] for r in rows if r.get("kind") == "enc"
                            and r.get("d") == 7 and r.get("bit") == 1 and r.get("K") == 32
                            and r.get("code", "bit") == "bit"), None)),
        ("e11 相位编码 d7 K32", "e11_memory",
         lambda rows: next((r["err"] for r in rows if r.get("kind") == "enc"
                            and r.get("d") == 7 and r.get("bit") == 1 and r.get("K") == 32
                            and r.get("code") == "phase"), None)),
        ("e10 最好星(后选)", "e10_code422",
         lambda rows: min((r["ps_err_perq"] for r in rows
                           if "ps_err_perq" in r and "kind" not in r), default=None)),
        ("e10 融合制备(FT验证)", "e10_code422",
         lambda rows: min((r["ps_err_perq"] for r in rows
                           if r.get("kind") == "fusion"), default=None)),
    ]
    print("\n多轮聚合(当日全部轮次,mean±CI / best / n):")
    for name, exp, ex in aggs:
        p = pool(exp, ex)
        if p.get("n"):
            print(f"  {name}: mean={p['mean']:.5f} CI[{p['lo']:.5f},{p['hi']:.5f}] "
                  f"best={p['best']:.5f} n={p['n']}")

    # e24 跨轮聚合(加强轮加密后 n=10/日):arms 在 summary 顶层,须读原始 dict;
    # fidelity/接受率越大越好 → best 取 pool_ci 的 worst 字段(=max)
    e24_docs = [json.loads(f.read_text())
                for f in sorted((R / "e24_logical_bell").glob("summary_Baihua_2*.json"))
                if TODAY in f.name]
    for name, ex24 in [
        ("e24 逻辑Bell F", lambda d: (d.get("arms", {}).get("logical") or {}).get("fidelity")),
        ("e24 物理Bell F", lambda d: (d.get("arms", {}).get("physical") or {}).get("fidelity")),
        ("e24 逻辑臂接受率", lambda d: (d.get("arms", {}).get("logical") or {}).get("acceptance")),
    ]:
        p = analysis.pool_ci([v for v in (ex24(d) for d in e24_docs)
                              if isinstance(v, (int, float))])
        if p.get("n"):
            print(f"  {name}: mean={p['mean']:.5f} CI[{p['lo']:.5f},{p['hi']:.5f}] "
                  f"best={p['worst']:.5f} n={p['n']}")


if __name__ == "__main__":
    main()
