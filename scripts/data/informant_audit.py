#!/usr/bin/env python3
"""Audit informant membership of the panel's multiz alignments against the
benchmark species panel.

For every leaf of each UCSC multiz tree, resolve the leaf to a scientific
name and NCBI taxid, then compare with ``benchmark/panel.tsv``:

* ``drop``      the leaf is a held-out panel species (same taxid): a training
                run must drop this row at cut time and declare it in
                ``alignment_rows_dropped`` (docs/benchmark.md section 3.3).
* ``train``     the leaf is a train panel species (the target or another
                train species); allowed.
* ``genus``     same genus as a panel species but a different species; kept,
                listed so the choice is visible.
* ``other``     no relation to the panel.
* ``unresolved`` the leaf could not be mapped to a taxid; listed for manual
                inspection.

Standard library only, Python 3.11. Nothing is downloaded by this script;
it reads the tree and name files fetched into ``--nh-dir`` beforehand. The
inputs, their URLs under https://hgdownload.soe.ucsc.edu/goldenPath/ and
https://api.genome.ucsc.edu/, and their checksums on 2026-09-18 are listed in
``informant_audit.inputs.md5`` next to this script; the file name in
``--nh-dir`` is the URL's basename, prefixed with ``<db>_`` where two tracks
would otherwise collide:

    mm39/multiz35way/mm39.35way.nh
    dm6/multiz124way/dm6.124way.sequenceNames.nh
    dm6/multiz124way/dm6.124way.scientificName.nh   -> dm6_dm6.124way.scientificName.nh
    dm6/multiz124way/dm6.124way.taxId.nh            -> dm6_dm6.124way.taxId.nh
    dm6/multiz27way/dm6.27way.nh
    ce11/multiz135way/ce11.135way.nh
    ce11/multiz135way/ce11.135way.scientificName.nh -> ce11_ce11.135way.scientificName.nh
    ce11/multiz135way/ce11.135way.taxId.nh          -> ce11_ce11.135way.taxId.nh
    sacCer3/multiz7way/7way.nh                      -> sacCer3_7way.nh
    hg38/multiz100way/hg38.100way.nh
    hg38/multiz100way/hg38.100way.scientificNames.nh
    galGal6/multiz77way/galGal6.77way.nh
    galGal6/multiz77way/galGal6.77way.scientificNames.nh
    https://api.genome.ucsc.edu/list/ucscGenomes    -> genomes.json

Resolution order per leaf: the track's own taxId/scientificName trees when
UCSC ships them (dm6 124-way, ce11 135-way); the UCSC genome list for active
database names; a small hand table for retired or non-UCSC leaves (casCan1,
GCF_003668045v3, the sacCer3 7-way sensu stricto yeasts); the same leaf name
resolved on another track; the track's scientificNames tree; and finally a
leaf whose name without its trailing version digits was resolved elsewhere.
The output lists the mechanism implicitly through taxid (empty when only the
name was resolved) and the relation column.

Result on 2026-09-18 (``scripts/data/informant_membership.tsv``): the
complete ``alignment_rows_dropped`` manifests for training are
``[hg38, galGal6]`` for mm39 35-way, ``[apiMel4]`` for dm6 124-way and
27-way, ``[ci3]`` for ce11 135-way and ``[]`` for sacCer3 7-way.
"""
from __future__ import annotations
import argparse, csv, json, re, sys
from pathlib import Path

LEAF = re.compile(r"(?:(?<=[(,\s])|^)([A-Za-z0-9_.\-]+)(?::-?[0-9.eE\-]+)?(?=\s*[,)])")

def leaves(newick: str) -> list[str]:
    s = re.sub(r"\s+", " ", newick.strip())
    # a leaf label is preceded by '(' or ',' (or whitespace inside a bare tree)
    out = []
    for m in re.finditer(r"([(,\s])\s*([A-Za-z0-9_.\-]+)", s):
        out.append(m.group(2))
    return out

def parallel(*files: Path) -> list[tuple[str, ...]]:
    cols = [leaves(f.read_text()) for f in files]
    assert len({len(c) for c in cols}) == 1, [len(c) for c in cols]
    return list(zip(*cols))

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nh-dir", required=True, help="directory with the fetched .nh / name files")
    ap.add_argument("--panel", default="benchmark/panel.tsv")
    ap.add_argument("--out", default="-")
    a = ap.parse_args()
    d = Path(a.nh_dir)
    genomes = json.loads((d / "genomes.json").read_text())["ucscGenomes"]
    db2 = {db: (g.get("scientificName", ""), int(g.get("taxId") or 0)) for db, g in genomes.items()}
    # hand-resolved leaves that are not UCSC databases
    db2.update({
        "GCF_003668045v3": ("Cricetulus griseus", 10029),   # mm39 35-way; NCBI assembly report GCF_003668045.3
        "casCan1": ("Castor canadensis", 51338),             # mm39 35-way; retired UCSC db, NCBI taxonomy 51338
        "sacPar": ("Saccharomyces paradoxus", 27291), "sacMik": ("Saccharomyces mikatae", 114525),
        "sacKud": ("Saccharomyces kudriavzevii", 114524), "sacBay": ("Saccharomyces bayanus", 4931),
        "sacCas": ("Naumovozyma castellii", 27288), "sacKlu": ("Lachancea kluyveri", 4934),
    })
    panel = list(csv.DictReader(open(a.panel), delimiter="\t"))
    by_tax = {int(r["tax_id"]): r for r in panel}
    by_genus = {}
    for r in panel:
        by_genus.setdefault(r["species"].split("_")[0], []).append(r)

    tracks = {
        "mm39/multiz35way":  ("train",   [("mm39.35way.nh", "db")]),
        "dm6/multiz124way":  ("train",   [("dm6.124way.sequenceNames.nh", "db"), ("dm6_dm6.124way.scientificName.nh", "sci"), ("dm6_dm6.124way.taxId.nh", "tax")]),
        "dm6/multiz27way":   ("train",   [("dm6.27way.nh", "db")]),
        "ce11/multiz135way": ("train",   [("ce11.135way.nh", "db"), ("ce11_ce11.135way.scientificName.nh", "sci"), ("ce11_ce11.135way.taxId.nh", "tax")]),
        "sacCer3/multiz7way": ("train",  [("sacCer3_7way.nh", "db")]),
        "hg38/multiz100way": ("eval",    [("hg38.100way.nh", "db"), ("hg38.100way.scientificNames.nh", "sci")]),
        "galGal6/multiz77way": ("eval",  [("galGal6.77way.nh", "db"), ("galGal6.77way.scientificNames.nh", "sci")]),
    }
    by_name = {r["species"].replace("_", " "): r for r in panel}
    resolved: dict[str, tuple[str, int]] = {}   # leaf -> (sci, taxid) from any track with a taxId file
    for track, (use, files) in tracks.items():
        kinds = [k for _, k in files]
        if "tax" in kinds:
            for rec in parallel(*[d / f for f, _ in files]):
                resolved[rec[kinds.index("db")]] = (rec[kinds.index("sci")].replace("_", " "), int(rec[kinds.index("tax")]))
    stem = lambda x: re.sub(r"\d+$", "", x)
    by_stem = {stem(k): v for k, v in resolved.items() if stem(k) != k}
    rows = []
    for track, (use, files) in tracks.items():
        recs = parallel(*[d / f for f, _ in files])
        kinds = [k for _, k in files]
        for rec in recs:
            leaf = rec[kinds.index("db")]
            sci, tax = "", 0
            if "tax" in kinds:
                tax = int(rec[kinds.index("tax")])
                sci = rec[kinds.index("sci")].replace("_", " ")
            elif leaf in db2:
                sci, tax = db2[leaf]
            elif leaf in resolved:
                sci, tax = resolved[leaf]
            elif "sci" in kinds:
                sci = rec[kinds.index("sci")].replace("_", " ")
            elif stem(leaf) in by_stem:
                sci, tax = by_stem[stem(leaf)]
            genus = sci.split(" ")[0] if sci else ""
            if tax in by_tax or (not tax and sci in by_name):
                p = by_tax[tax] if tax else by_name[sci]
                rel = "drop" if p["split"].startswith("heldout") else "train"
                panel_sp, split = p["species"], p["split"]
            elif genus and genus in by_genus:
                rel = "genus"; panel_sp = ";".join(r["species"] for r in by_genus[genus]); split = ";".join(r["split"] for r in by_genus[genus])
            elif tax or sci:
                rel, panel_sp, split = "other", "", ""
            else:
                rel, panel_sp, split = "unresolved", "", ""
            rows.append({"track": track, "use": use, "leaf": leaf, "scientific_name": sci, "taxid": tax or "",
                         "relation": rel, "panel_species": panel_sp, "panel_split": split})
    out = sys.stdout if a.out == "-" else open(a.out, "w")
    w = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rows)
    # summary to stderr
    from collections import Counter
    for track in tracks:
        c = Counter(r["relation"] for r in rows if r["track"] == track)
        n = sum(c.values())
        flagged = [f'{r["leaf"]}={r["panel_species"]}' for r in rows if r["track"] == track and r["relation"] in ("drop", "train", "genus")]
        print(f"{track}: {n} leaves; " + ", ".join(f"{k}={v}" for k, v in sorted(c.items())) + (" | " + "; ".join(flagged) if flagged else ""), file=sys.stderr)
    return 0

if __name__ == "__main__":
    sys.exit(main())
