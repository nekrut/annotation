#!/usr/bin/env python3
"""Merge the five Phase 1 review bibliographies into one deduplicated file.

Input:  relay/artifacts/T-human-{002,003,004,005,012}/refs.bib
Output: docs/refs/refs.bib

Deduplication key is the normalized DOI. Where several reviews carry the same
DOI, fields are merged field-by-field, preferring the longest non-empty value,
and the entry records which reviews contributed it. Entries without a DOI are
kept only if they carry a URL, and are keyed by normalized title.

Standard library only, Python 3.11. Run from the repository root:

    python3 scripts/review/merge_refs.py
"""

from __future__ import annotations

import csv
import io
import re
import sys
import unicodedata
from pathlib import Path

CHECK_DATE = "2026-09-10"

SOURCES = {
    "T-human-002": "lenin",
    "T-human-003": "engels",
    "T-human-004": "stalin",
    "T-human-005": "marx",
    "T-human-012": "trotsky",
}

FIELD_ORDER = [
    "author",
    "title",
    "journal",
    "booktitle",
    "howpublished",
    "year",
    "volume",
    "number",
    "pages",
    "doi",
    "url",
    "pmid",
    "pmcid",
    "eprint",
    "publisher",
    "note",
    "verified",
]

# Institutional or corporate authors, where "last word of the author string"
# makes a poor citation key.
KEY_OVERRIDES = {
    "information2026egapx": "ncbi2026egapx",
}

STOPWORDS = {
    "a", "an", "and", "the", "of", "for", "in", "on", "with", "to", "by",
    "from", "using", "via", "new", "at",
}


def strip_braces(value: str) -> str:
    """Remove wrapping braces, but only when the first brace is the partner of
    the last one. `{BRAKER3}: ... {TSEBRA}` must survive intact."""
    value = value.strip()
    while len(value) > 1 and value.startswith("{") and value.endswith("}"):
        depth = 0
        for i, ch in enumerate(value):
            if ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0:
                    break
        if i != len(value) - 1:
            break
        value = value[1:-1].strip()
    return value


def parse_fields(body: str) -> dict:
    """Read `name = value` pairs from an entry body.

    Values are scanned with a brace counter rather than a regex: author lists
    in these bibliographies nest braces several levels deep (`{{Br{\\r{u}}na}}`)
    and a fixed-depth pattern silently truncates them.
    """
    fields = {}
    pos = 0
    while True:
        m = re.compile(r"(\w+)\s*=\s*").search(body, pos)
        if not m:
            break
        name = m.group(1).lower()
        i = m.end()
        if i < len(body) and body[i] == "{":
            depth = 0
            j = i
            while j < len(body):
                if body[j] == "{":
                    depth += 1
                elif body[j] == "}":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            raw, pos = body[i:j + 1], j + 1
        elif i < len(body) and body[i] == '"':
            j = body.find('"', i + 1)
            j = len(body) - 1 if j == -1 else j
            raw, pos = body[i + 1:j], j + 1
        else:
            j = body.find(",", i)
            j = len(body) if j == -1 else j
            raw, pos = body[i:j], j + 1
        fields[name] = " ".join(strip_braces(raw.strip().rstrip(",")).split())
    return fields


def parse_bib(text: str) -> list[dict]:
    """A tolerant BibTeX reader: entries, keys and brace/quote-delimited fields."""
    entries = []
    i = 0
    while True:
        at = text.find("@", i)
        if at == -1:
            break
        m = re.match(r"@(\w+)\s*\{\s*([^,]+),", text[at:])
        if not m:
            i = at + 1
            continue
        entry_type, key = m.group(1).lower(), m.group(2).strip()
        # find the matching closing brace of the entry, starting at the brace
        # that opens it (the one right after the entry type)
        depth = 0
        start_body = text.index("{", at)
        j = start_body
        while j < len(text):
            if text[j] == "{":
                depth += 1
            elif text[j] == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        body = text[start_body + 1:j]
        entries.append(
            {"type": entry_type, "key": key, "fields": parse_fields(body)}
        )
        i = j + 1
    return entries


def norm_doi(doi: str) -> str:
    doi = doi.strip().lower()
    doi = re.sub(r"^(https?://)?(dx\.)?doi\.org/", "", doi)
    doi = re.sub(r"^doi:\s*", "", doi)
    return doi.rstrip(".")


def find_url(fields: dict) -> str:
    """A URL for entries with no DOI: an explicit url, or one embedded in
    howpublished / note (several reviews record repository documentation that
    way)."""
    if fields.get("url"):
        return fields["url"].strip()
    for name in ("howpublished", "note"):
        m = re.search(r"https?://\S+", fields.get(name, ""))
        if m:
            return m.group(0).rstrip("}.,")
    return ""


def norm_url(url: str) -> str:
    """Collapse a GitHub URL to owner/repo so that a README, a pinned blob and
    the repository root are recognized as one source."""
    url = url.strip().rstrip("/")
    m = re.match(r"https?://(?:www\.)?github\.com/([^/]+)/([^/]+)", url)
    if m:
        return "github.com/{}/{}".format(m.group(1), m.group(2))
    return re.sub(r"^https?://(?:www\.)?", "", url).lower()


def norm_title(title: str) -> str:
    title = unicodedata.normalize("NFKD", title.lower())
    return re.sub(r"[^a-z0-9]+", "", title)


def first_surname(author: str) -> str:
    if not author:
        return "anon"
    first = author.split(" and ")[0].strip()
    if "," in first:
        surname = first.split(",")[0]
    else:
        surname = first.split()[-1] if first.split() else "anon"
    surname = unicodedata.normalize("NFKD", surname)
    return re.sub(r"[^a-z]", "", surname.lower()) or "anon"


def title_word(title: str) -> str:
    for word in re.findall(r"[A-Za-z]+", title):
        low = word.lower()
        if low not in STOPWORDS and len(low) > 2:
            return low
    return "untitled"


def make_key(fields: dict) -> str:
    return "{}{}{}".format(
        first_surname(fields.get("author", "")),
        fields.get("year", "nd"),
        title_word(fields.get("title", "")),
    )


def merge_fields(a: dict, b: dict) -> dict:
    out = dict(a)
    for name, value in b.items():
        if name == "note":
            continue
        if len(value) > len(out.get(name, "")):
            out[name] = value
    return out


def format_entry(key: str, entry: dict) -> str:
    fields = entry["fields"]
    ordered = [f for f in FIELD_ORDER if fields.get(f)]
    ordered += sorted(f for f in fields if f not in FIELD_ORDER and fields[f])
    width = max(len(f) for f in ordered)
    lines = ["@{}{{{},".format(entry["type"], key)]
    for name in ordered:
        lines.append("  {:<{w}} = {{{}}},".format(name, fields[name], w=width))
    lines.append("}")
    return "\n".join(lines)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    merged: dict[str, dict] = {}
    counts = {}
    dropped = []
    for task, agent in SOURCES.items():
        path = root / "relay" / "artifacts" / task / "refs.bib"
        entries = parse_bib(path.read_text(encoding="utf-8"))
        counts[agent] = len(entries)
        for entry in entries:
            fields = entry["fields"]
            doi = norm_doi(fields.get("doi", ""))
            url = find_url(fields)
            if doi:
                ident = "doi:" + doi
                fields["doi"] = doi
            elif url:
                ident = "url:" + norm_url(url)
                fields.setdefault("url", url)
            else:
                dropped.append((agent, entry["key"]))
                continue
            if ident in merged:
                merged[ident]["fields"] = merge_fields(
                    merged[ident]["fields"], fields
                )
                merged[ident]["sources"].add(agent)
            else:
                merged[ident] = {
                    "type": entry["type"],
                    "fields": fields,
                    "sources": {agent},
                }

    # assign collision-free citation keys
    keyed: dict[str, dict] = {}
    for ident, entry in merged.items():
        key = make_key(entry["fields"])
        key = KEY_OVERRIDES.get(key, key)
        if key in keyed:
            suffix = "b"
            while key + suffix in keyed:
                suffix = chr(ord(suffix) + 1)
            key = key + suffix
        keyed[key] = entry

    # If a verification pass has already run, carry its verdict into each
    # entry so that a DOI Crossref resolves to the wrong paper cannot be cited
    # from this file by accident.
    verdicts = {}
    verification = root / "docs" / "review" / "doi-verification.tsv"
    if verification.exists():
        text = verification.read_text(encoding="utf-8")
        records = list(csv.reader(io.StringIO(text), delimiter="\t"))
        header, body = records[0], records[1:]
        width = len(header)
        bad = [i + 2 for i, r in enumerate(body) if len(r) != width]
        if bad:
            # A record split across physical lines would silently hand a
            # truncated title, or an agent name, to the DO NOT CITE notice
            # below. Refuse the file instead of writing a wrong verdict.
            raise SystemExit(
                f"{verification.name}: {len(bad)} record(s) are not {width} fields "
                f"wide (lines {bad[:5]}); re-run scripts/review/verify_dois.py"
            )
        for cols in body:
            verdicts[cols[0]] = (cols[2], cols[3])

    out_lines = [
        "% Merged Phase 1 bibliography for the eukaryotic gene prediction charter.",
        "% Produced by scripts/review/merge_refs.py from the five independent",
        "% review bibliographies under relay/artifacts/. Do not edit by hand:",
        "% edit the source bibliography or the merge script and re-run.",
        "%",
        "% Source entry counts: "
        + ", ".join(f"{a} {counts[a]}" for a in sorted(counts)),
        f"% Merged unique works: {len(keyed)}",
        "% Deduplication key: normalized DOI (fallback: normalized title + URL).",
        "% The `reviews` field lists which independent reviews carried the work.",
        "",
    ]
    for key in sorted(keyed):
        entry = keyed[key]
        entry["fields"]["reviews"] = ", ".join(sorted(entry["sources"]))
        status, registered = verdicts.get(key, ("", ""))
        if status == "mismatch":
            entry["fields"]["verified"] = (
                "DO NOT CITE: Crossref resolves this DOI to a different work, "
                f"\"{registered}\" (checked {CHECK_DATE})"
            )
        elif status == "unresolved":
            entry["fields"]["verified"] = (
                f"DO NOT CITE: Crossref has no record of this DOI (checked {CHECK_DATE})"
            )
        elif status:
            entry["fields"]["verified"] = f"crossref {status}, {CHECK_DATE}"
        out_lines.append(format_entry(key, entry))
        out_lines.append("")

    out_path = root / "docs" / "refs" / "refs.bib"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(out_lines), encoding="utf-8")

    with_doi = sum(1 for e in keyed.values() if e["fields"].get("doi"))
    shared = sum(1 for e in keyed.values() if len(e["sources"]) > 1)
    print(f"read: {counts}")
    print(f"wrote {len(keyed)} entries to {out_path.relative_to(root)}")
    print(f"  with DOI: {with_doi}; found by >1 review: {shared}")
    if dropped:
        print(f"  dropped (no DOI and no URL): {dropped}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
