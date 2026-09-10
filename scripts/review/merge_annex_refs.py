#!/usr/bin/env python3
"""Collect the works cited only in the `engels` and `stalin` annex messages.

Both agents submitted a base review plus 22 annex messages, and both covers
state the annexes take precedence over the base draft. The annexes cite works
their own `refs.bib` never got: `engels`'s cover says "eight distinct DOI
sources there are absent from the base bibliography" and `stalin`'s index
lists twelve. Those two counts are self-reported and overlap; this script
measures the union instead.

It scans every message from `engels` or `stalin` carrying `task: T-human-003`
or `task: T-human-004`, extracts DOIs, subtracts the DOIs already present in
the five artifact bibliographies, fetches Crossref metadata for what is left,
and writes:

    docs/review/annex-dois.tsv   one row per annex-only DOI, with provenance
    docs/review/annex-refs.bib   BibTeX for the resolvable ones

`merge_refs.py` reads `annex-refs.bib` as a sixth source, so the annex works
land in `docs/refs/refs.bib` with `reviews = {engels (annex)}` and an `annex`
field naming the messages that cite them.

Crossref is the only network access, and only with --fetch. Without it the
script re-reads the existing annex-refs.bib and only rewrites the TSV, so the
diff is reproducible offline.

Standard library only, Python 3.11. Run from the repository root:

    python3 scripts/review/merge_annex_refs.py [--fetch]
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
from merge_refs import SOURCES, format_entry, norm_doi, parse_bib  # noqa: E402

UA = (
    "relay-review-synthesis/0.1 "
    "(https://github.com/nekrut/annotation; mailto:anton@nekrut.org)"
)
CHECK_DATE = "2026-09-10"
ANNEX_TASKS = {"T-human-003": "engels", "T-human-004": "stalin"}
COLUMNS = ("doi", "agents", "annexes", "status", "title", "first_annex")

# A DOI in running prose picks up whatever punctuation follows it. Stop at
# whitespace, a closing bracket, a backtick, or sentence punctuation.
DOI_RE = re.compile(r"10\.\d{4,9}/[^\s\"<>()\[\]{}`,;]+", re.I)


def clean_doi(raw: str) -> str:
    """Trim trailing prose punctuation and the version/format suffixes that
    bioRxiv links carry (`...434297v3`, `...05449-z.pdf`)."""
    doi = norm_doi(raw).rstrip(".,;:")
    doi = re.sub(r"\.(pdf|full|supplementary)$", "", doi)
    if re.match(r"10\.1101/2\d{3}\.", doi):
        doi = re.sub(r"v\d+$", "", doi)
    return doi


def annex_messages(root: Path) -> list[tuple[str, str, str]]:
    """(agent, message id, text) for every engels/stalin annex message."""
    out = []
    for path in sorted((root / "relay" / "messages").glob("*.md")):
        text = path.read_text(encoding="utf-8")
        sender = re.search(r"^from: (\S+)\s*$", text, re.M)
        task = re.search(r"^task: (\S+)\s*$", text, re.M)
        if not sender or not task:
            continue
        agent, task_id = sender.group(1), task.group(1)
        if ANNEX_TASKS.get(task_id) != agent:
            continue
        out.append((agent, path.stem, text))
    return out


def known_dois(root: Path) -> set[str]:
    """Every DOI the five artifact bibliographies already carry."""
    known = set()
    for task in SOURCES:
        path = root / "relay" / "artifacts" / task / "refs.bib"
        for entry in parse_bib(path.read_text(encoding="utf-8")):
            doi = norm_doi(entry["fields"].get("doi", ""))
            if doi:
                known.add(doi)
    return known


def crossref(doi: str) -> dict | None:
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="/")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=30) as fh:
            return json.load(fh)["message"]
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return None
        raise


def flat(value: str) -> str:
    """Crossref titles carry embedded newlines; a TSV record must be one line."""
    return " ".join(str(value).split())


def bib_from_crossref(msg: dict, doi: str) -> dict:
    authors = []
    for person in msg.get("author", []) or []:
        if person.get("family"):
            given = person.get("given", "").strip()
            authors.append(
                f"{person['family']}, {given}" if given else person["family"]
            )
        elif person.get("name"):
            authors.append("{" + person["name"] + "}")
    container = (msg.get("container-title") or [""])[0]
    issued = msg.get("issued", {}).get("date-parts", [[None]])[0]
    fields = {
        "author": " and ".join(authors),
        "title": flat((msg.get("title") or [""])[0]),
        "journal": flat(container),
        "year": str(issued[0]) if issued and issued[0] else "",
        "volume": flat(msg.get("volume", "")),
        "number": flat(msg.get("issue", "")),
        "pages": flat(msg.get("page", "")),
        "doi": doi,
        "url": msg.get("URL", ""),
    }
    entry_type = "article" if container else "misc"
    if msg.get("type") == "posted-content":
        fields["journal"] = fields["journal"] or "bioRxiv"
        entry_type = "article"
    return {"type": entry_type, "fields": {k: v for k, v in fields.items() if v}}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--fetch",
        action="store_true",
        help="query Crossref for DOIs not already in annex-refs.bib",
    )
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    known = known_dois(root)

    cited: dict[str, dict[str, set[str]]] = {}
    for agent, msg_id, text in annex_messages(root):
        for raw in DOI_RE.findall(text):
            doi = clean_doi(raw)
            if not doi or doi in known:
                continue
            record = cited.setdefault(doi, {"agents": set(), "annexes": set()})
            record["agents"].add(agent)
            record["annexes"].add(msg_id)

    out_bib = root / "docs" / "review" / "annex-refs.bib"
    existing = {}
    if out_bib.exists():
        for entry in parse_bib(out_bib.read_text(encoding="utf-8")):
            existing[norm_doi(entry["fields"].get("doi", ""))] = entry

    resolved: dict[str, dict] = {}
    status: dict[str, str] = {}
    for doi in sorted(cited):
        if doi in existing:
            resolved[doi] = existing[doi]
            status[doi] = "cached"
            continue
        if not args.fetch:
            status[doi] = "not-fetched"
            continue
        msg = crossref(doi)
        time.sleep(1.0)
        if msg is None:
            status[doi] = "not-in-crossref"
            continue
        resolved[doi] = bib_from_crossref(msg, doi)
        status[doi] = "ok"

    lines = [
        "% Works cited only in the engels and stalin annex messages.",
        "% Produced by scripts/review/merge_annex_refs.py --fetch; metadata is",
        f"% Crossref's, retrieved {CHECK_DATE}. Do not edit by hand.",
        "%",
        f"% Annex messages scanned: {len(annex_messages(root))}",
        f"% DOIs cited in annexes but absent from the five artifact",
        f"% bibliographies: {len(cited)}; resolvable: {len(resolved)}",
        "",
    ]
    for doi in sorted(resolved):
        entry = resolved[doi]
        record = cited[doi]
        entry["fields"]["annex"] = "{}: {}".format(
            ", ".join(sorted(record["agents"])),
            ", ".join(sorted(record["annexes"])),
        )
        entry["fields"]["reviews"] = ", ".join(
            f"{a} (annex)" for a in sorted(record["agents"])
        )
        key = entry.get("key") or ""
        if not key:
            from merge_refs import make_key

            key = make_key(entry["fields"])
        entry["key"] = key
        lines.append(format_entry(key, entry))
        lines.append("")
    out_bib.write_text("\n".join(lines), encoding="utf-8")

    out_tsv = root / "docs" / "review" / "annex-dois.tsv"
    buf = io.StringIO()
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
    writer.writerow(COLUMNS)
    for doi in sorted(cited):
        entry = resolved.get(doi)
        writer.writerow(
            [
                doi,
                ",".join(sorted(cited[doi]["agents"])),
                len(cited[doi]["annexes"]),
                status.get(doi, ""),
                flat(entry["fields"]["title"]) if entry else "",
                sorted(cited[doi]["annexes"])[0],
            ]
        )
    text = buf.getvalue()
    out_tsv.write_text(text, encoding="utf-8")
    bad = [
        i + 2
        for i, row in enumerate(list(csv.reader(io.StringIO(text), delimiter="\t"))[1:])
        if len(row) != len(COLUMNS)
    ]
    if bad:
        raise SystemExit(f"{out_tsv.name}: records not {len(COLUMNS)} wide: {bad[:5]}")

    print(f"annex-only DOIs: {len(cited)}; resolved: {len(resolved)}")
    for doi in sorted(cited):
        if status.get(doi) not in ("ok", "cached"):
            print(f"  {status.get(doi)}: {doi}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
