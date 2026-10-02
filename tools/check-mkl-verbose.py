import re
import sys
from pathlib import Path

# Each DGEMM is attributed to its requested thread budget. Aggregate counts
# alone can hide a call exceeding its own limit or a serial public pipeline.
seen = {}
current = None
for line in Path(sys.argv[1]).read_text(encoding="utf-8").splitlines():
    label = re.match(r"^(THREAD_PROBE|PUBLIC_PIPELINE)\s+([12])\s*$", line)
    if label:
        current = (label.group(1), int(label.group(2)))
        seen.setdefault(current, [])
    elif current and "DGEMM" in line and "MKL_VERBOSE" in line:
        match = re.search(r"NThr:\s*(\d+)", line)
        if not match:
            raise SystemExit("DGEMM verbose line is missing NThr")
        threads = int(match.group(1))
        if threads < 1 or threads > current[1]:
            raise SystemExit(f"DGEMM exceeded its call budget: {current}, NThr={threads}")
        seen[current].append(threads)
for kind in ("THREAD_PROBE", "PUBLIC_PIPELINE"):
    for budget in (1, 2):
        values = seen.get((kind, budget), [])
        if not values or budget not in values:
            raise SystemExit(f"Missing DGEMM evidence for {kind} with NThr={budget}")
print("Verified oneMKL per-call DGEMM teams:", seen)
