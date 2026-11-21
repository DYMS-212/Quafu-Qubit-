"""等待某台后端从 Calibrating/维护恢复(当天周期重试)。

平台标定/维护可能中途修好;一次探到不可用就整天放弃太浪费。本工具在当天
窗口内周期重探 client.status():
  - 一旦可用(队列深度为 int)→ 打印放行,退出码 0
  - 到封顶 max_min 仍不可用 → 打印放弃(交明天 cron),退出码 1
  - max_min=0 → 只探一次(供"整道跑完后判断是否需重试补跑")

用法: ./qrun tools/wait_available.py <chip> [max_min=300] [poll_min=20]
自检: ./qrun tools/wait_available.py --selftest   # 纯离线,不打真机
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import client  # noqa: E402


def wait_available(chip: str, max_min: float = 300.0, poll_min: float = 20.0) -> int:
    """周期探测直到 chip 可用(→0)或达封顶仍不可用(→1)。max_min=0 只探一次。"""
    deadline = time.time() + max_min * 60
    attempt = 0
    while True:
        attempt += 1
        try:
            v = client.status().get(chip)
        except Exception as e:  # noqa: BLE001
            v = f"探测异常:{e}"
        if isinstance(v, int) and not isinstance(v, bool):   # 队列深度=可用
            print(f"[{time.strftime('%H:%M')}] {chip} 可用(队列={v}),第 {attempt} 次探测→放行")
            return 0
        left = deadline - time.time()
        if left <= 0:
            print(f"[{time.strftime('%H:%M')}] {chip}={v},探 {attempt} 次达封顶 "
                  f"{max_min:.0f} 分钟仍未恢复→放弃(交明天 cron)")
            return 1
        print(f"[{time.strftime('%H:%M')}] {chip}={v}(平台标定/维护中),{poll_min:.0f} "
              f"分钟后重试(第 {attempt} 次,剩 {left / 60:.0f} 分钟窗口)")
        time.sleep(min(poll_min * 60, left))


def _selftest() -> int:
    """离线自检:mock client.status,校验 max_min=0 对不可用即返 1、可用返 0。"""
    orig = client.status
    try:
        client.status = lambda: {"X": "Calibrating"}
        assert wait_available("X", max_min=0, poll_min=1) == 1, "不可用应返回 1"
        client.status = lambda: {"X": True}     # bool 不算可用队列深度
        assert wait_available("X", max_min=0, poll_min=1) == 1, "bool 不应判为可用"
        client.status = lambda: {"X": 7}
        assert wait_available("X", max_min=0, poll_min=1) == 0, "可用应返回 0"
        print("wait_available selftest OK")
        return 0
    finally:
        client.status = orig


def main() -> int:
    if len(sys.argv) > 1 and sys.argv[1] == "--selftest":
        return _selftest()
    chip = sys.argv[1] if len(sys.argv) > 1 else "Shenglian"
    max_min = float(sys.argv[2]) if len(sys.argv) > 2 else 300.0
    poll_min = float(sys.argv[3]) if len(sys.argv) > 3 else 20.0
    return wait_available(chip, max_min, poll_min)


if __name__ == "__main__":
    sys.exit(main())
