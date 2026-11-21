"""补收:扫描 results/ 里"已提交但没收到结果"的记录,按 tid 重新拉取并落盘。

用途:批次被打断(保险丝/超时/进程退出)后,任务仍在平台队列里跑完了,
     用本脚本把结果找回来。安全:只读平台结果,不新提交任何任务。

运行: ./qrun tools/recollect.py [天数]   # 默认只看最近 2 天的记录
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import client  # noqa: E402

PROJ = Path(__file__).resolve().parent.parent
DAYS = int(sys.argv[1]) if len(sys.argv) > 1 else 2


def main():
    cutoff = time.time() - DAYS * 86400
    missing = []
    for f in (PROJ / "results").rglob("*__*.json"):
        if f.stat().st_mtime < cutoff:
            continue
        try:
            r = json.loads(f.read_text())
        except Exception:  # noqa: BLE001
            continue
        tid = r.get("tid")
        if isinstance(tid, int) and not (r.get("result") or {}).get("count"):
            missing.append((f, r))
    print(f"扫描完成:{len(missing)} 条缺结果记录(近 {DAYS} 天)")
    ok = fail = 0
    for f, r in missing:
        try:
            res = client.result(r["tid"])
            if isinstance(res, dict) and res.get("count"):
                r["result"] = res
                r["fetched_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
                f.write_text(json.dumps(r, ensure_ascii=False, indent=2))
                ok += 1
                print(f"  ✓ {r.get('label')} ({f.parent.name})")
            else:
                st = res.get("status") if isinstance(res, dict) else res
                fail += 1
                print(f"  ✗ {r.get('label')}: {st}")
        except Exception as e:  # noqa: BLE001
            fail += 1
            print(f"  ✗ {r.get('label')}: {str(e)[:80]}")
    print(f"补收 {ok} 条,仍缺 {fail} 条(Failed 的任务不会有结果,属正常)")


if __name__ == "__main__":
    main()
