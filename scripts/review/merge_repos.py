#!/usr/bin/env python3
"""Merge the five Phase 1 repository inventories into one table.

The five reviews used different column names and, in two cases, different
methods for counting commits (GitHub REST `since` filters versus a shallow
`git clone --shallow-since`). All five snapshots were taken on 2026-09-09, so
where their numbers differ the difference is method, not drift, and that is
worth keeping visible rather than averaging away.

Output: docs/review/repos.tsv, one row per repository, with the range of
values each field took across the reviews that recorded it and a
`disagreements` column naming the fields where they did not agree.

Standard library only. Run from the repository root:

    python3 scripts/review/merge_repos.py
"""

from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from merge_refs import SOURCES  # noqa: E402

# canonical name -> the column names the five reviews used for it
COLUMNS = {
    "last_commit": ["last_commit", "last_commit_date"],
    "commits_12m": ["commits_12m", "commits_last_12_months", "commits_last_12mo",
                    "commits_last_12m"],
    "open_issues": ["open_issues"],
    "stars": ["stars"],
    "language": ["language"],
    "licence": ["licence", "license"],
    "install": ["install_tested", "install_test", "readme_install_fresh_env",
                "readme_install_worked"],
}
NUMERIC = {"commits_12m", "open_issues", "stars"}
# stars and issues move on their own; only flag a numeric gap this large
NUMERIC_TOLERANCE = 0.10


def slug(url: str) -> str:
    url = url.strip().rstrip("/")
    m = re.search(r"github\.com/([^/]+)/([^/\s]+)", url)
    if m:
        return "{}/{}".format(m.group(1), m.group(2))
    return url


def pick(row: dict, names: list[str]) -> str:
    for name in names:
        value = (row.get(name) or "").strip()
        if value and value not in {"?", "-", "n/a", "NA"}:
            return value
    return ""


def numeric_disagreement(values: list[str]) -> bool:
    nums = []
    for value in values:
        m = re.search(r"\d+", value.replace(",", ""))
        if m:
            nums.append(int(m.group(0)))
    if len(nums) < 2:
        return False
    lo, hi = min(nums), max(nums)
    return hi - lo > max(1, NUMERIC_TOLERANCE * hi)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    repos: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for task, agent in SOURCES.items():
        path = root / "relay" / "artifacts" / task / "repos.tsv"
        with path.open(encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                url = (row.get("url") or "").strip()
                if not url:
                    continue
                repos[slug(url)][agent] = {
                    field: pick(row, names) for field, names in COLUMNS.items()
                }

    # Two annex messages carry a repository snapshot that never reached their
    # author's own repos.tsv: `abacus-gene/paml` (stalin 0016) and
    # `Jstacs/Jstacs`, which hosts GeMoMa (engels 0018). They are transcribed
    # into docs/review/annex-repos.tsv with the annex id in the `annex` column,
    # and enter the merge as "<agent> (annex)".
    annex_path = root / "docs" / "review" / "annex-repos.tsv"
    annex_rows = 0
    if annex_path.exists():
        with annex_path.open(encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                url = (row.get("url") or "").strip()
                if not url:
                    continue
                label = "{} (annex)".format((row.get("agent") or "annex").strip())
                repos[slug(url)][label] = {
                    field: pick(row, names) for field, names in COLUMNS.items()
                }
                annex_rows += 1

    out_rows = []
    for repo in sorted(repos, key=str.lower):
        seen = repos[repo]
        row = {"repo": repo, "reviews": ",".join(sorted(seen))}
        disagree = []
        for field in COLUMNS:
            values = [v[field] for v in seen.values() if v[field]]
            uniq = sorted(set(values))
            row[field] = " | ".join(uniq) if uniq else ""
            if len(uniq) > 1:
                if field in NUMERIC:
                    if numeric_disagreement(values):
                        disagree.append(field)
                elif field in {"licence", "language"}:
                    disagree.append(field)
        row["disagreements"] = ",".join(disagree)
        out_rows.append(row)

    out = root / "docs" / "review" / "repos.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    header = ["repo", "reviews", *COLUMNS, "disagreements"]
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(out_rows)

    multi = sum(1 for r in out_rows if len(r["reviews"].split(",")) > 1)
    flagged = [r["repo"] for r in out_rows if r["disagreements"]]
    print(f"wrote {len(out_rows)} repositories to {out.relative_to(root)}")
    print(f"  annex snapshots folded in: {annex_rows}")
    print(f"  seen by more than one review: {multi}")
    print(f"  fields disagreeing across reviews: {len(flagged)} repos")
    for repo in flagged:
        row = next(r for r in out_rows if r["repo"] == repo)
        print(f"  ! {repo}: {row['disagreements']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
