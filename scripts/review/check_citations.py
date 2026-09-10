#!/usr/bin/env python3
"""Check that every citation key used in docs/review/ exists in docs/refs/refs.bib,
and report bibliography entries nothing cites.

Citation form in these documents is `[key]` or `[key1; key2]`, where a key is
surname + year + word, as produced by merge_refs.py.

    python3 scripts/review/check_citations.py

Exits 1 if a cited key is missing.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from merge_refs import parse_bib  # noqa: E402

KEY = r"[a-z][a-z0-9]+\d{4}[a-z]+"
CITE = re.compile(rf"\[({KEY}(?:;\s*{KEY})*)\]")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    bib = root / "docs" / "refs" / "refs.bib"
    keys = {e["key"] for e in parse_bib(bib.read_text(encoding="utf-8"))}

    cited: dict[str, set[str]] = {}
    for path in sorted((root / "docs" / "review").glob("*.md")):
        for m in CITE.finditer(path.read_text(encoding="utf-8")):
            for key in re.split(r";\s*", m.group(1)):
                cited.setdefault(key, set()).add(path.name)

    missing = sorted(k for k in cited if k not in keys)
    unused = sorted(keys - set(cited))
    print(f"{len(keys)} entries in {bib.relative_to(root)}; "
          f"{len(cited)} distinct keys cited; {len(unused)} uncited")
    for key in missing:
        print(f"  MISSING {key} (cited in {', '.join(sorted(cited[key]))})")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
