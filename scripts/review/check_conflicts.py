#!/usr/bin/env python3
"""Report works that two or more reviews cite under different DOIs.

Runs over the five Phase 1 bibliographies and groups entries whose titles are
the same work: identical after normalization, or sharing at least
JACCARD_MIN of their content words. A group carrying more than one distinct
DOI is a citation conflict: the same work has been recorded with different
identifiers, so at most one of them can be right.

The output is the evidence behind docs/review/disagreements.md section 1.
Run from the repository root:

    python3 scripts/review/check_conflicts.py
"""

from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from merge_refs import SOURCES, STOPWORDS, norm_doi, parse_bib  # noqa: E402

JACCARD_MIN = 0.6


def content_words(title: str) -> frozenset[str]:
    words = re.findall(r"[a-z0-9]+", title.lower())
    return frozenset(w for w in words if w not in STOPWORDS and len(w) > 2)


def group_titles(rows: list[tuple]) -> list[list[tuple]]:
    """Single-link clustering of (agent, doi, title) rows by title overlap."""
    groups: list[list[tuple]] = []
    keys: list[set[frozenset[str]]] = []
    for row in rows:
        words = content_words(row[2])
        target = None
        for i, group_keys in enumerate(keys):
            for other in group_keys:
                union = words | other
                if union and len(words & other) / len(union) >= JACCARD_MIN:
                    target = i
                    break
            if target is not None:
                break
        if target is None:
            groups.append([row])
            keys.append({words})
        else:
            groups[target].append(row)
            keys[target].add(words)
    return groups


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    all_rows = []
    for task, agent in SOURCES.items():
        path = root / "relay" / "artifacts" / task / "refs.bib"
        for entry in parse_bib(path.read_text(encoding="utf-8")):
            title = entry["fields"].get("title", "")
            doi = norm_doi(entry["fields"].get("doi", ""))
            if title and doi:
                all_rows.append((agent, doi, " ".join(title.split())))

    conflicts = 0
    for rows in group_titles(all_rows):
        dois = {doi for _, doi, _ in rows}
        if len(dois) < 2:
            continue
        conflicts += 1
        print(f"\n{rows[0][2]}")
        tally = defaultdict(list)
        for agent, doi, _ in rows:
            tally[doi].append(agent)
        for doi in sorted(tally, key=lambda d: (-len(tally[d]), d)):
            print(f"  {doi:<40} {', '.join(sorted(tally[doi]))}")

    print(f"\n{conflicts} title(s) cited under more than one DOI")
    return 0


if __name__ == "__main__":
    sys.exit(main())
