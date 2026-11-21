"""汇总层身份+溯源戳:每个 summary 顶层嵌一份,让跨日/跨轮统计可归因、可复现、可审计。

背景:原始 record(client._save_record)有物理比特/tid,但各实验写 summary 时普遍
只留聚合指标、把身份和溯源丢了 → 跨日 p_L(d) 误差棒会把"不同物理链、不同标定日"
的点当同一条盲池("串味"),且锚点同名原地覆盖后无法复现当天选了谁。stamp() 统一补齐。
"""
from __future__ import annotations

import time


def stamp(anchor: dict | None = None, qubits=None, tids=None,
          chip: str = "Baihua", extra: dict | None = None) -> dict:
    """构造 summary 顶层 provenance 段。
    anchor: anchors.latest()/load() 返回的锚点 dict(含 generated_from/date/
            chip_calibration_time/score_mode);headline 选址实验必传。
    qubits: 本实验实际用的物理比特(链/星/面板);None 则由各数据行自带。
    tids:   本实验提交的任务 id 列表(审计凭据)。
    """
    a = anchor or {}
    p = {
        "chip": chip,
        "anchor_source": a.get("generated_from"),        # selfcal_eta / vendor_readout / composite
        "anchor_score_mode": a.get("score_mode"),        # readout / gate / composite
        "anchor_date": a.get("date"),                    # 当日锚点日期;冻结锚点为 None
        "chip_calibration_time": a.get("chip_calibration_time"),
        "run_time": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    if qubits is not None:
        p["qubits"] = qubits
    if tids is not None:
        p["tids"] = list(tids)
    if extra:
        p.update(extra)
    return p


def tids_of(records) -> list:
    """从一批 submit_job/run_batch 返回的 record 里抽 tid(int),供 stamp(tids=...)。"""
    out = []
    for r in records or []:
        t = r.get("tid") if isinstance(r, dict) else None
        if isinstance(t, int):
            out.append(t)
    return out
