"""与 Quafu 真机交互。

已验证事实(勿改):
  - token 在工作区根 ../.env,key = QUAFU_TOKEN;Task 优先读环境变量 QPU_API_TOKEN
  - Task() 构造即联网 verify,平台偶发抖动 → 构造必须重试
  - shots 必须是 1024 的整数倍;每天 1000 任务
  - count 的 bit 串序 = c_{n-1}..c1 c0(最右 = c0)
  - target_qubits[i] = 逻辑 q_i 映射到的物理比特;只能放 available 比特,否则 transpiler 拒
  - 并发跑两个批次易触发 ConnectionReset → 实验串行单批
"""
from __future__ import annotations

import json
import os
import threading
import time
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
# 工作区(放 .env)= 项目上一级;可用 QFAB_WORKSPACE 覆盖(跨机移植:Mac=夸父/, 服务器=kuafu/)
WORKSPACE = Path(os.environ.get("QFAB_WORKSPACE", PROJ.parent))
RESULTS = PROJ / "results"

DEFAULT_OPTIONS = {
    "compiler": "quarkcircuit",
    "correct": False,
    "open_dd": None,
    "target_qubits": [],
}

_FINISHED = {"Finished", "Failed", "Completed", "Cancelled", "Done"}


# ----------------------------------------------------------------- 基础设施
def _call_timeout(fn, *a, _fuse=60.0, **kw):
    """熔断:在 daemon 线程里跑 fn,超 _fuse 秒未返回则抛 TimeoutError。
    quark SDK 的 status/result/run/Task() 都无 timeout,平台一 hang 就永久阻塞
    (曾把主道锁死 9.5h)。此熔断把"卡住不返回"转成可捕获的异常 → _retry 重试
    → 最终仍超时则抛,由上层 ||true 跳过该实验,主道不再被单个调用锁死。
    超时后底层 SDK 线程仍 hang 但为 daemon,进程退出时自动清理(轻微线程泄漏,无害)。"""
    box: dict = {}

    def _run():
        try:
            box["ok"] = fn(*a, **kw)
        except BaseException as e:  # noqa: BLE001
            box["err"] = e

    t = threading.Thread(target=_run, daemon=True)
    t.start()
    t.join(_fuse)
    if t.is_alive():
        raise TimeoutError(f"SDK 调用超 {_fuse:.0f}s 未返回(平台 hang),熔断跳过")
    if "err" in box:
        raise box["err"]
    return box.get("ok")


def _retry(fn, *a, tries=6, delay=2.0, **kw):
    last = None
    for i in range(tries):
        try:
            return fn(*a, **kw)
        except Exception as e:  # noqa: BLE001
            last = e
            if i < tries - 1:
                time.sleep(delay * (i + 1))
    raise last


# 多 token:QFAB_TOKEN_ID=A(默认→QUAFU_TOKEN)/B(→QUAFU_TOKEN2)/C(→QUAFU_TOKEN3)...
# 每个 token 平台各有 1000/天额度 → 多 key=多额度;额度按 token 分别计数(见 _quota_file)。
TOKEN_ID = os.environ.get("QFAB_TOKEN_ID", "A").strip().upper()


def _token_key() -> str:
    n = ord(TOKEN_ID) - ord("A")  # A=0→QUAFU_TOKEN, B=1→QUAFU_TOKEN2, ...
    return "QUAFU_TOKEN" if n == 0 else f"QUAFU_TOKEN{n + 1}"


def load_token() -> str:
    key_name = _token_key()
    alt = ("QPU_API_TOKEN", "token") if TOKEN_ID == "A" else ()
    raw = (WORKSPACE / ".env").read_text(encoding="utf-8", errors="ignore")
    for line in raw.splitlines():
        key, _, val = line.strip().partition("=")
        if key.strip() == key_name or key.strip() in alt:
            tok = val.strip().strip('"').strip("'")
            os.environ["QPU_API_TOKEN"] = tok
            return tok
    raise RuntimeError(f"在 .env 里没找到 {key_name}(QFAB_TOKEN_ID={TOKEN_ID})")


_task = None


def get_task():
    global _task
    if _task is None:
        from quark import Task
        _task = _retry(lambda: _call_timeout(Task, load_token(), _fuse=60),
                       tries=8, delay=4)
    return _task


def status(tid: int = 0):
    """tid=0 时返回各后端排队状态(安全、不耗额度)。熔断 60s 防平台 hang 锁死主道。"""
    return _retry(lambda: _call_timeout(get_task().status, tid, _fuse=60), tries=4)


def result(tid: int, timeout: float = 0.0):
    # 熔断 45s:collect 里高频轮询,单个 tid 查询 hang 不再卡死整个 collect,
    # deadline 逻辑才能真正生效(否则 while 循环卡在这行,轮不到检查超时)。
    return _retry(lambda: _call_timeout(get_task().result, tid, timeout, _fuse=45),
                  tries=3)


# ----------------------------------------------------------------- 芯片信息
_chip_cache: dict = {}

# 网页控制台同源端点(无鉴权 GET,schema 与 SDK Backend.chip_info 逐字段一致,已核对)。
# quarkcircuit 0.5.5 的 Backend 有芯片白名单(Baihua/Dongling/Miaofeng),Shenglian 被拒
# → 白名单外机器走此端点;拓扑/标定/priority_qubits/坐标全量可得。
CHIP_DETAILS_URL = "https://quafu-sqc.baqis.ac.cn/sqc-api/chip-details/{chip}"


def _chip_info_http(chip: str) -> dict:
    import requests
    r = requests.get(CHIP_DETAILS_URL.format(chip=chip), timeout=30)
    r.raise_for_status()
    info = r.json()
    if not isinstance(info, dict) or "qubits_info" not in info:
        raise RuntimeError(f"chip-details 返回异常: {str(info)[:200]}")
    return info


def get_chip_info(chip: str = "Baihua", refresh: bool = False) -> dict:
    """拉取并缓存芯片标定(不耗额度),同时按日期落盘一份供漂移追踪。"""
    if not refresh and chip in _chip_cache:
        return _chip_cache[chip]
    try:
        get_task()  # 确保 token 已注入环境
        from quark.circuit import Backend
        info = Backend(chip).chip_info
        if not isinstance(info, dict) or "qubits_info" not in info:
            raise RuntimeError("SDK chip_info 为空")
    except Exception:
        info = _retry(_chip_info_http, chip, tries=4, delay=3)
    _chip_cache[chip] = info
    d = RESULTS / "chip_snapshots"
    d.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y-%m-%d")
    (d / f"chip_{chip}_{stamp}.json").write_text(
        json.dumps(info, ensure_ascii=False, indent=2), encoding="utf-8")
    return info


# ----------------------------------------------------------------- 提交/收取
def _save_record(record: dict) -> Path:
    d = RESULTS / record["exp"]
    d.mkdir(parents=True, exist_ok=True)
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in record["label"])
    path = d / f"{safe}__{record.get('tid', 'noid')}.json"
    path.write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
    return path


# 无人为额度帽子:真实每日上限由平台决定(隐藏,不一定是 1000)。RUNAWAY_GUARD 仅防脚本 bug 无限提交。
RUNAWAY_GUARD = int(os.environ.get("QFAB_RUNAWAY_GUARD", "100000"))


def _quota_file(day: str | None = None) -> Path:
    day = day or time.strftime("%Y-%m-%d")
    prefix = "" if TOKEN_ID == "A" else f"{TOKEN_ID}_"  # token A 保持旧命名(向后兼容)
    return RESULTS / "quota_log" / f"{prefix}{day}.txt"


def quota_used(day: str | None = None) -> int:
    """当前 token 当日"成功被平台接收"的任务数(不计被拒/失败的提交)。"""
    f = _quota_file(day)
    return sum(1 for _ in f.open()) if f.exists() else 0


def quota_used_all(day: str | None = None) -> int:
    """所有 token 当日成功提交总数(A+B+...)。"""
    day = day or time.strftime("%Y-%m-%d")
    d = RESULTS / "quota_log"
    tot = 0
    for tid in _available_token_ids():
        f = d / (f"{day}.txt" if tid == "A" else f"{tid}_{day}.txt")
        if f.exists():
            tot += sum(1 for _ in f.open())
    return tot


def _quota_log_success():
    (RESULTS / "quota_log").mkdir(parents=True, exist_ok=True)
    with _quota_file().open("a") as fh:
        fh.write(time.strftime("%H:%M:%S") + "\n")


# ---- 多 token failover:当前 token 被平台拒(额度耗尽)或故障 → 自动切下一个 ----
_TOKEN_IDS: list[str] | None = None
_ALL_EXHAUSTED = False


def _available_token_ids() -> list[str]:
    """.env 里配了哪些 token:A=QUAFU_TOKEN, B=QUAFU_TOKEN2, C=QUAFU_TOKEN3 ..."""
    global _TOKEN_IDS
    if _TOKEN_IDS is None:
        raw = (WORKSPACE / ".env").read_text(encoding="utf-8", errors="ignore")
        keys = {ln.split("=", 1)[0].strip() for ln in raw.splitlines()}
        _TOKEN_IDS = [tid for i, tid in enumerate("ABCDEFGH")
                      if ("QUAFU_TOKEN" if i == 0 else f"QUAFU_TOKEN{i + 1}") in keys] or ["A"]
    return _TOKEN_IDS


def _failover() -> bool:
    """切到下一个可用 token;成功 True,已是最后一个 False。"""
    global TOKEN_ID, _task
    ids = _available_token_ids()
    if TOKEN_ID in ids and ids.index(TOKEN_ID) + 1 < len(ids):
        nxt = ids[ids.index(TOKEN_ID) + 1]
        print(f"⚠️ token {TOKEN_ID} 提交失败(额度耗尽/故障)→ 自动切 token {nxt}")
        TOKEN_ID = nxt
        _task = None  # 强制用新 token 重新鉴权
        return True
    return False


def _submit_with_failover(task: dict):
    """提交 task;被平台拒/故障则自动切下一个 token 重试。返回 tid(int)或 error(dict)。"""
    global _ALL_EXHAUSTED
    if _ALL_EXHAUSTED:
        return {"error": "all tokens exhausted today"}
    if quota_used() >= RUNAWAY_GUARD:  # 纯防脚本跑飞,非人为额度帽
        raise RuntimeError(f"RUNAWAY_GUARD {RUNAWAY_GUARD} 触发(疑似 bug)")
    while True:
        try:
            tid = _retry(lambda: _call_timeout(get_task().run, task, _fuse=90),
                         tries=5, delay=3)
        except Exception as e:  # noqa: BLE001
            tid = {"error": repr(e)}
        if not isinstance(tid, dict):
            _quota_log_success()  # 只记真正被接收的
            return tid
        if _failover():           # 当前 token 用尽/故障 → 切下一个重试
            continue
        _ALL_EXHAUSTED = True      # 所有 token 都用尽,本进程后续直接短路
        return tid


def submit_job(circuit: str, chip: str = "Baihua", shots: int = 2048,
               name: str = "qfab", options: dict | None = None,
               exp: str = "misc", label: str | None = None,
               meta: dict | None = None) -> dict:
    """异步提交,立即返回 record(result=None)并落盘。无人为额度帽:
    跑到平台自己拒 → 自动 failover 到下一个 token(A→B→…),全用尽才停。"""
    assert shots % 1024 == 0, "shots 必须是 1024 的整数倍"
    opts = dict(DEFAULT_OPTIONS)
    if options:
        opts.update(options)
    task = {"chip": chip, "name": name, "circuit": circuit,
            "shots": shots, "options": opts}
    tid = _submit_with_failover(task)
    record = {
        "exp": exp, "label": label or name, "chip": chip, "shots": shots,
        "tid": tid, "token": TOKEN_ID, "submitted_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "task": task, "meta": meta or {}, "result": None,
    }
    if isinstance(tid, dict):  # 所有 token 都提交失败
        record["error"] = tid
    _save_record(record)
    return record


def _is_done(res) -> bool:
    return isinstance(res, dict) and (res.get("status") in _FINISHED or bool(res.get("count")))


def collect(records: list[dict], timeout: float = 1200, poll: float = 3.0,
            verbose: bool = True) -> list[dict]:
    """轮询一批 record,完成即落盘。"""
    pending = {r["tid"]: r for r in records
               if r.get("tid") is not None and not isinstance(r["tid"], dict)}
    total, done = len(pending), 0
    deadline = time.time() + timeout
    while pending and time.time() < deadline:
        for tid in list(pending):
            try:
                res = result(tid)
            except Exception:  # noqa: BLE001
                continue
            if _is_done(res):
                r = pending.pop(tid)
                r["result"] = res
                r["fetched_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
                _save_record(r)
                done += 1
                if verbose:
                    print(f"  [{done}/{total}] {r['label']} -> {res.get('status')}")
        if pending:
            time.sleep(poll)
    for tid, r in pending.items():
        try:
            r["result"] = result(tid)
        except Exception:  # noqa: BLE001
            pass
        _save_record(r)
        if verbose:
            print(f"  [超时未完成] {r['label']} tid={tid}")
    return records


def run_batch(specs: list[dict], timeout: float = 1200, poll: float = 3.0) -> list[dict]:
    """specs 每项 = submit_job 的关键字参数。一次性异步提交 → 统一轮询。"""
    print(f"批量提交 {len(specs)} 个任务 ...")
    recs = [submit_job(**s) for s in specs]
    bad = [r for r in recs if isinstance(r.get("tid"), dict)]
    if bad:
        for r in bad:
            print(f"  [提交失败] {r['label']}: {r['tid']}")
    print(f"提交成功 {len(recs) - len(bad)}/{len(specs)},开始轮询 ...")
    return collect(recs, timeout=timeout, poll=poll)


def counts(record: dict) -> dict:
    r = record.get("result") or {}
    return r.get("count", {}) if isinstance(r, dict) else {}
