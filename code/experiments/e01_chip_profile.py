"""e01 · 芯片画像:拓扑/可用性/厂家发布标定概览(不耗额度)。

运行: ./qrun experiments/e01_chip_profile.py [chip]
"""
from __future__ import annotations

import json
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qfab import chip as chipmod  # noqa: E402

CHIP = sys.argv[1] if len(sys.argv) > 1 else "Baihua"
TODAY = time.strftime("%Y-%m-%d")
OUT = Path(__file__).resolve().parent.parent / "results" / "e01_chip"


def main():
    ch = chipmod.load(CHIP, refresh=True)
    av = ch.available_qubits()
    deg = Counter(len([x for x in ch.neighbors(q) if x in av]) for q in av)
    stars = [q for q in av
             if len([x for x in ch.neighbors(q)
                     if x in av and (ch.edge_fid(q, x) or 0) > 0.5]) >= 3]
    profile = {
        "date": TODAY, "chip": CHIP, **ch.summary(),
        "available": av,
        "degree_distribution": dict(sorted(deg.items())),
        "n_degree3_stars": len(stars),
        "unavailable_qubits": sorted(set(ch.qubits) - set(av)),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"profile_{CHIP}_{TODAY}.json").write_text(
        json.dumps(profile, ensure_ascii=False, indent=2))
    print(json.dumps(profile, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
