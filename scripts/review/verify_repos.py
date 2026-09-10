#!/usr/bin/env python3
"""Re-check the repository rows the five reviews disagreed about.

docs/review/repos.tsv flags repositories where the reviews recorded different
licences, languages or commit counts, and the merge cannot tell which is
right. This asks GitHub directly for the repositories in that flagged set,
plus any slug passed on the command line, and writes the answer to
docs/review/repo-verification.tsv.

The unauthenticated GitHub API allows 60 requests an hour, which is why this
checks only the disputed rows rather than all of them. `gh` is used when it is
installed and logged in, which lifts the limit.

    python3 scripts/review/verify_repos.py [extra/slug ...]
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

UA = "relay-review-synthesis/0.1 (mailto:anton@nekrut.org)"
FIELDS = ["full_name", "stargazers_count", "open_issues_count", "language",
          "archived", "pushed_at"]


def fetch(slug: str) -> tuple[str, dict]:
    if shutil.which("gh"):
        proc = subprocess.run(["gh", "api", f"repos/{slug}"],
                              capture_output=True, text=True)
        if proc.returncode == 0:
            return "ok", json.loads(proc.stdout)
        if "Not Found" in proc.stderr:
            return "404", {}
    req = urllib.request.Request("https://api.github.com/repos/" + slug,
                                 headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=20) as fh:
            # urllib follows the redirect GitHub issues for a renamed
            # repository, so a differing full_name below means "renamed".
            return "ok", json.load(fh)
    except urllib.error.HTTPError as exc:
        return str(exc.code), {}
    except Exception as exc:
        return type(exc).__name__, {}


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parents[2]
    merged = root / "docs" / "review" / "repos.tsv"
    slugs = list(argv)
    if merged.exists():
        with merged.open(encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["disagreements"] and row["repo"] not in slugs:
                    slugs.append(row["repo"])

    rows = []
    for slug in slugs:
        status, data = fetch(slug)
        licence = ((data.get("license") or {}).get("spdx_id") or "-") if data else "-"
        rows.append({
            "queried": slug,
            "status": status,
            "resolves_to": data.get("full_name", "-"),
            "renamed": "yes" if data and data.get("full_name", "").lower() != slug.lower() else "no",
            "stars": data.get("stargazers_count", "-"),
            "open_issues": data.get("open_issues_count", "-"),
            "language": data.get("language") or "-",
            "spdx_licence": licence,
            "archived": data.get("archived", "-"),
            "pushed_at": data.get("pushed_at", "-"),
        })
        time.sleep(0.3)

    out = root / "docs" / "review" / "repo-verification.tsv"
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"checked {len(rows)} repositories -> {out.relative_to(root)}")
    for row in rows:
        if row["status"] != "ok" or row["renamed"] == "yes":
            print(f"  ! {row['queried']} status={row['status']} "
                  f"resolves_to={row['resolves_to']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
