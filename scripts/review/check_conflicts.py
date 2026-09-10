#!/usr/bin/env python3
"""Report works that two or more reviews cite under different DOIs.

Runs over the five Phase 1 bibliographies and groups entries whose titles are
the same work: identical after normalization, or sharing at least
JACCARD_MIN of their content words. A group carrying more than one distinct DOI is reported, but not every such
group is an error. Two categories come out separately:

  conflict     the same work recorded under two identifiers that cannot both
               be right; this is the evidence behind section 1 of
               docs/review/disagreements.md
  preprint_of  one preprint DOI (bioRxiv, arXiv, Research Square) and one
               journal DOI for the same work -- both identifiers are valid
               and refer to different versions. Nucleotide Transformer
               (10.1101/2023.01.11.523679 / 10.1038/s41592-024-02523-z) and
               SegmentNT (10.1101/2024.03.14.584712 / 10.1038/s41592-025-02881-2)
               are the two in this bibliography.

Only the conflict count feeds section 1. Reading the whole list as errors
attributes a legitimate preprint/journal pair to whichever review recorded
the preprint.
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

# bioRxiv/medRxiv share the 10.1101 prefix with Cold Spring Harbor journals
# (Genome Research is 10.1101/gr.NNNNNN), so the date-stamped path, not the
# prefix, is what identifies a preprint.
PREPRINT = (
    re.compile(r"^10\.1101/\d{4}\.\d{2}\.\d{2}\."),   # bioRxiv, medRxiv
    re.compile(r"^10\.48550/"),                          # arXiv
    re.compile(r"^10\.21203/"),                          # Research Square
)


def is_preprint(doi: str) -> bool:
    return any(rx.match(doi) for rx in PREPRINT)


def classify(dois: set[str]) -> str:
    """A preprint plus exactly one journal DOI is two versions, not a clash."""
    published = {d for d in dois if not is_preprint(d)}
    if len(published) <= 1 and len(dois) - len(published) >= 1:
        return "preprint_of"
    return "conflict"


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

    found = {"conflict": [], "preprint_of": []}
    for rows in group_titles(all_rows):
        dois = {doi for _, doi, _ in rows}
        if len(dois) < 2:
            continue
        found[classify(dois)].append(rows)

    for kind in ("conflict", "preprint_of"):
        if not found[kind]:
            continue
        print(f"\n== {kind} ({len(found[kind])})")
        for rows in found[kind]:
            print(f"\n{rows[0][2]}")
            tally = defaultdict(list)
            for agent, doi, _ in rows:
                tally[doi].append(agent)
            for doi in sorted(tally, key=lambda d: (-len(tally[d]), d)):
                print(f"  {doi:<40} {', '.join(sorted(tally[doi]))}")

    print(f"\n{len(found['conflict'])} title(s) cited under conflicting DOIs; "
          f"{len(found['preprint_of'])} preprint/journal pair(s), which are not errors")
    return 0


if __name__ == "__main__":
    sys.exit(main())
