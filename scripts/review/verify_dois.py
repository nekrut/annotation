#!/usr/bin/env python3
"""Resolve every DOI in the merged bibliography against Crossref.

For each entry in docs/refs/refs.bib this fetches api.crossref.org/works/<doi>
and compares the registered title with the title the reviews recorded. The
result is written to docs/review/doi-verification.tsv with one row per entry:

    key  doi  status  registered_title  title_overlap  reviews

status is `ok` (resolves, titles agree), `mismatch` (resolves to a different
work), `unresolved` (Crossref has no record for a DOI it should serve),
`not-in-crossref` (a prefix Crossref does not serve at all: arXiv DOIs are
registered with DataCite), `no-title` (Crossref has the record but no title
to compare, common for journal supplements), or `no-doi`.

Every field is whitespace-flattened before it is written. Crossref titles
carry embedded newlines (the registered Helixer title breaks after "Gene"),
and an embedded newline splits one logical record across several physical
TSV lines, which is what any reader -- `csv.reader`, `cut`, `awk` -- sees.
The file is read back after writing and the run fails if any record is not
exactly six fields wide or if the key column does not match the bibliography.

Standard library only. Run from the repository root:

    python3 scripts/review/verify_dois.py [--limit N]
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_conflicts import content_words  # noqa: E402
from merge_refs import parse_bib  # noqa: E402

# Crossref asks for a contact address in the User-Agent so it can reach the
# caller about excessive traffic; anton is this runner's operator.
UA = "relay-review-synthesis/0.1 (https://github.com/nekrut/annotation; mailto:anton@nekrut.org)"
OVERLAP_OK = 0.6
COLUMNS = ("key", "doi", "status", "registered_title", "title_overlap", "reviews")


def flatten(value: object) -> str:
    """One logical record per physical line: no tabs, no newlines, no runs."""
    return re.sub(r"\s+", " ", str(value)).strip()


def serialize(rows: list[tuple]) -> str:
    out = ["\t".join(COLUMNS)]
    for row in rows:
        if len(row) != len(COLUMNS):
            raise ValueError(f"row has {len(row)} fields, expected {len(COLUMNS)}: {row!r}")
        out.append("\t".join(flatten(v) for v in row))
    return "\n".join(out) + "\n"


def check_serialization(text: str, rows: list[tuple]) -> None:
    """Parse the written file the way a consumer does and fail on any drift."""
    records = list(csv.reader(io.StringIO(text), delimiter="\t"))
    if len(records) != len(rows) + 1:
        raise ValueError(f"{len(records)} parsed records for {len(rows)} entries plus header")
    bad = [i for i, r in enumerate(records) if len(r) != len(COLUMNS)]
    if bad:
        raise ValueError(f"records with wrong width at lines {[i + 1 for i in bad]}")
    if tuple(records[0]) != COLUMNS:
        raise ValueError(f"header is {records[0]}")
    written = [r[0] for r in records[1:]]
    if written != [flatten(r[0]) for r in rows]:
        raise ValueError("key column does not match the bibliography order")


def strip_markup(title: str) -> str:
    return re.sub(r"<[^>]+>", "", title)


def crossref(doi: str) -> tuple[str, str]:
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="/")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=20) as fh:
            message = json.load(fh)["message"]
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            # arXiv registers with DataCite, so Crossref legitimately has no
            # record; that is not evidence against the DOI.
            return ("not-in-crossref" if doi.startswith("10.48550/") else "unresolved", "")
        return (f"http-{exc.code}", "")
    except Exception as exc:  # network failure, timeout, bad payload
        return ("error", type(exc).__name__)
    titles = message.get("title") or []
    return "found", strip_markup(titles[0]) if titles else ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="stop after N entries")
    ap.add_argument("--delay", type=float, default=0.4, help="seconds between calls")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    entries = parse_bib((root / "docs" / "refs" / "refs.bib").read_text(encoding="utf-8"))
    if args.limit:
        entries = entries[: args.limit]

    rows = []
    for entry in entries:
        key = entry["key"]
        fields = entry["fields"]
        doi = fields.get("doi", "")
        reviews = fields.get("reviews", "")
        recorded = fields.get("title", "")
        if not doi:
            rows.append((key, "-", "no-doi", "", "", reviews))
            continue
        state, registered = crossref(doi)
        time.sleep(args.delay)
        if state != "found":
            rows.append((key, doi, state if state != "error" else "error",
                         registered, "", reviews))
            continue
        if not registered:
            rows.append((key, doi, "no-title", "", "", reviews))
            continue
        a, b = content_words(recorded), content_words(registered)
        overlap = len(a & b) / len(a | b) if (a | b) else 0.0
        status = "ok" if overlap >= OVERLAP_OK else "mismatch"
        rows.append((key, doi, status, registered, f"{overlap:.2f}", reviews))

    out = root / "docs" / "review" / "doi-verification.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    text = serialize(rows)
    check_serialization(text, rows)
    out.write_text(text, encoding="utf-8")

    tally = {}
    for row in rows:
        tally[row[2]] = tally.get(row[2], 0) + 1
    print(f"checked {len(rows)} entries -> {out.relative_to(root)}")
    for status in sorted(tally):
        print(f"  {status}: {tally[status]}")
    for row in rows:
        if row[2] in ("mismatch", "unresolved"):
            print(f"  ! {row[0]} {row[1]} {row[2]} {row[3][:60]} [{row[5]}]")
    return 0


if __name__ == "__main__":

    sys.exit(main())
