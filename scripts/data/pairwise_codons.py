#!/usr/bin/env python3
"""Export reference-versus-informant codon pairs from a fetched window, for
the KA/KS comparative baseline (T-human-010) and any other pairwise codon
model.

Standard library only (Python 3.11).  Reads the ``<stem>.*`` files written by
``fetch_window.py`` and reuses the MAF, tree and isoform code of
``cut_windows.py``.

Why a separate export exists
----------------------------
codeml's pairwise mode (``runmode = -2``) forces ``cleandata = 1`` and its
reader deletes every codon column in which *any* input sequence carries a
gap or an ambiguity, before the pair is fitted (relay note
20260909T164106Z-stalin-0016, from a trace of ``abacus-gene/paml`` at
``4c7902f``).  Handing codeml a multispecies alignment therefore removes
codons that are clean in the target/informant pair.  This script
materialises each pair on its own, from the reference CDS and one informant
row, and records for every informant how many codons the CDS has, how many
survive the *pairwise* rule, how many would survive the *complete-case* rule
over all rows, and why each dropped codon was dropped.  ``--complete-case``
switches the written pairs to the complete-case columns so the two policies
can be fitted on identical inputs and compared.

Codon rules
-----------
A reference codon is three CDS bases in transcription order (reverse
complemented on the minus strand; a codon may span a splice junction).  The
frame is taken from the CDS start, which is exact for complete RefSeq and
Ensembl transcripts; ``cds_length_mod3`` in the sidecar is non-zero for a
partial CDS and the trailing partial codon is dropped.  A reference terminal
stop codon is removed (``reference_stop_removed``), because codeml rejects
stops.  A stop inside the reference CDS (selenocysteine, readthrough, or a
mis-set frame) is likewise excluded and its index listed in
``reference_internal_stops``; ``codons`` counts what remains and
``codons_in_cds`` what the CDS had.  For one informant a codon is *retained* when all three aligned bases
are A, C, G or T with no informant insertion between two consecutive codon
bases.  Otherwise it is dropped with the first applicable reason:
``unaligned`` (no MAF row covers the base), ``gap``, ``ambiguous`` (N or
another letter), ``insertion``, or ``informant_stop`` (the informant codon is
TAA, TAG or TGA; kept with ``--keep-stops``).  Identity is the fraction of
retained codons in which the two codons are equal.

Outputs (``<out>/<name>.<transcript>.*``)::

  .<informant>.phy   sequential PHYLIP, two rows, retained codons only
  .pairs.tsv         one row per informant: distance, codon counts by fate,
                     identity, and the complete-case count
  .sidecar.json      parameters, transcript, checksums of the inputs, the
                     per-codon fate matrix summary and the same table as
                     the TSV
  .all.phy           with --write-multi: every retained row over the
                     complete-case columns (the input that would reproduce
                     codeml's own deletion), for the comparison

Examples::

    python3 scripts/data/pairwise_codons.py --stem /tmp/win/Adh/Adh_124 \
        --out /tmp/pairs --isoforms longest-cds --drop-species apiMel4
    python3 scripts/data/pairwise_codons.py --stem /tmp/win/Adh/Adh_124 \
        --out /tmp/pairs --transcript NM_057262.3 --informants droSim1,droYak1
    python3 scripts/data/pairwise_codons.py --self-test
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cut_windows as cw  # noqa: E402

TOOL_VERSION = "0.1"
STOPS = {"TAA", "TAG", "TGA"}
REASONS = ("unaligned", "gap", "ambiguous", "insertion", "informant_stop")
CODE = "ACGT"


def load_window(stem: str):
    """Manifest, reference name, window, sequence, MAF blocks, trees."""
    with open(stem + ".manifest.json") as fh:
        manifest = json.load(fh)
    ref = manifest.get("assembly") or manifest.get("species")
    if manifest.get("source") == "ensembl":
        ref = manifest.get("species", ref)
    ref = ref.split(".")[0]
    win = manifest["window"]
    start, end = int(win["start"]), int(win["end"])
    with open(stem + ".fa") as fh:
        seq = "".join(ln.strip() for ln in fh if not ln.startswith(">")).upper()
    if len(seq) != end - start:
        raise SystemExit(f"{stem}: FASTA length {len(seq)} != window {end - start}")
    with open(stem + ".maf") as fh:
        blocks = cw.maf_blocks(fh.read())
    trees = []
    if os.path.exists(stem + ".nh"):
        with open(stem + ".nh") as fh:
            trees = [t.strip() + ";" for t in fh.read().split(";") if t.strip()]
    return manifest, ref, win["chrom"], start, end, seq, blocks, trees


def informant_rows(manifest: dict, ref: str, blocks, trees) -> tuple[list[str], dict[str, float]]:
    """Informant names in the cutter's row order (ancestral rows excluded)
    and their patristic distance from the reference."""
    sources: list[str] = []
    for b in blocks:
        for r in b:
            nm = cw.maf_source(r[1])
            if nm not in sources:
                sources.append(nm)
    declared = set(manifest.get("ancestral_sequences") or [])

    def is_ancestor(nm: str) -> bool:
        return nm in declared or cw.ANCESTOR_RE.match(nm) is not None or nm.lower().startswith("ancestor")

    names = set(sources) | {ref}
    dist: dict[str, float] = {}
    leaves_by_tree = []
    for t in trees:
        lv, d = cw.newick_leaves_and_distances(t, ref, names)
        leaves_by_tree.append(lv)
        for k, v in d.items():
            dist.setdefault(k, v)
    members = manifest.get("species_set_members")
    if members:
        leaves = sorted(members)
    elif len(trees) == 1:
        leaves = leaves_by_tree[0]
    else:
        leaves = sorted(s for s in sources if not is_ancestor(s))
    leaves = [n for n in leaves if n != ref and n.split(".")[0] != ref and not is_ancestor(n)]
    others = [s for s in sources if s not in leaves and s != ref and not is_ancestor(s)]
    return leaves + others, dist


def cds_codons(t: dict, start: int, end: int) -> tuple[list[tuple[int, int, int]], int, list[int]]:
    """Genomic positions of each codon's bases in transcription order,
    ``cds_length % 3`` and the CDS positions outside the window."""
    pos: list[int] = []
    for a, b in sorted(t.get("cds", [])):
        pos.extend(range(a, b))
    if t.get("strand") == "-":
        pos.reverse()
    outside = [q for q in pos if q < start or q >= end]
    codons = [(pos[i], pos[i + 1], pos[i + 2]) for i in range(0, len(pos) - len(pos) % 3, 3)]
    return codons, len(pos) % 3, outside


def ref_codon(seq: str, start: int, codon, strand: str) -> str:
    s = "".join(seq[q - start] for q in codon)
    return s.translate(cw.COMP) if strand == "-" else s


def informant_codon(inf: bytearray, ins: bytearray, start: int, codon, strand: str) -> tuple[str, str]:
    """(codon string, reason) with reason '' when retained by the pairwise rule."""
    vals = [inf[q - start] for q in codon]
    if any(v == cw.INF_UNALIGNED for v in vals):
        return "", "unaligned"
    if any(v == cw.INF_GAP for v in vals):
        return "", "gap"
    if any(v == cw.INF_OTHER for v in vals):
        return "", "ambiguous"
    lo = sorted(codon)
    for g in (lo[0], lo[1]):
        if g + 1 in lo and ins[g - start] > 0:
            return "", "insertion"
    s = "".join(CODE[v] for v in vals)
    if strand == "-":
        s = s.translate(cw.COMP)
    return s, ""


def phylip(rows: list[tuple[str, str]]) -> str:
    n = len(rows[0][1])
    out = [f"{len(rows)} {n}"]
    for name, s in rows:
        out.append(f"{name}  {s}")
    return "\n".join(out) + "\n"


def export(stem: str, out_dir: str, transcript: str | None, isoforms: str, informants: list[str] | None,
           drop: set[str], complete_case: bool, keep_stops: bool, write_multi: bool,
           transcript_types: str = "benchmark", quiet: bool = False) -> list[dict]:
    manifest, ref, chrom, start, end, seq, blocks, trees = load_window(stem)
    name = os.path.basename(stem)
    with open(stem + ".annotation.json") as fh:
        ann = json.load(fh)
    transcripts = ann.get("transcripts", ann if isinstance(ann, list) else [])
    transcripts, _ = cw.select_transcripts(transcripts, transcript_types)
    if transcript:
        chosen = [t for t in transcripts if str(t.get("id")) == transcript
                  or cw._strip_version(str(t.get("id"))) == cw._strip_version(transcript)]
        if not chosen:
            raise SystemExit(f"transcript {transcript} not in {stem}.annotation.json")
    else:
        chosen, _, _ = cw.select_isoforms(transcripts, isoforms)
    chosen = [t for t in chosen if t.get("cds")]
    rows, dist = informant_rows(manifest, ref, blocks, trees)
    rows = [r for r in rows if r not in drop]
    if informants:
        missing = [r for r in informants if r not in rows]
        if missing:
            raise SystemExit(f"informants not in the window's rows: {missing}")
        rows = [r for r in rows if r in informants]
    inf, ins, block_stats = cw.paint_informants(blocks, ref, start, end, rows)
    checks = {}
    for ext in (".maf", ".fa", ".annotation.json", ".manifest.json"):
        with open(stem + ext, "rb") as fh:
            checks[ext] = hashlib.sha256(fh.read()).hexdigest()
    os.makedirs(out_dir, exist_ok=True)
    results = []
    for t in chosen:
        tid = str(t.get("id"))
        strand = t.get("strand", "+")
        codons, mod3, outside = cds_codons(t, start, end)
        if outside:
            print(f"{name}/{tid}: {len(outside)} CDS bases outside the window; skipped", file=sys.stderr)
            continue
        refc = [ref_codon(seq, start, c, strand) for c in codons]
        stop_removed = False
        if refc and refc[-1] in STOPS:
            codons, refc = codons[:-1], refc[:-1]
            stop_removed = True
        codons_in_cds = len(codons)
        internal_stops = [i for i, c in enumerate(refc) if c in STOPS]
        if internal_stops and not keep_stops:
            keep = [i for i in range(len(codons)) if i not in set(internal_stops)]
            codons, refc = [codons[i] for i in keep], [refc[i] for i in keep]
        # fate[k][i]: ('' retained | reason), codon string
        fates: list[list[tuple[str, str]]] = []
        for k in range(len(rows)):
            row = []
            for c in codons:
                s, why = informant_codon(inf[k], ins[k], start, c, strand)
                if not why and s in STOPS and not keep_stops:
                    s, why = "", "informant_stop"
                row.append((s, why))
            fates.append(row)
        n = len(codons)
        cc_cols = [i for i in range(n) if all(not fates[k][i][1] for k in range(len(rows)))]
        cc_set = set(cc_cols)
        table = []
        for k, r in enumerate(rows):
            kept = [i for i in range(n) if not fates[k][i][1]]
            counts = {why: sum(1 for i in range(n) if fates[k][i][1] == why) for why in REASONS}
            use = cc_cols if complete_case else kept
            ident = sum(1 for i in use if fates[k][i][0] == refc[i])
            entry = {"informant": r, "dist": dist.get(r), "codons": n,
                     "retained_pairwise": len(kept), "retained_complete_case": len(cc_cols),
                     "written": len(use), "identical": ident,
                     "identity": (ident / len(use)) if use else None, **counts}
            table.append(entry)
            if use:
                rs = "".join(refc[i] for i in use)
                qs = "".join(fates[k][i][0] for i in use)
                with open(os.path.join(out_dir, f"{name}.{tid}.{r}.phy"), "w") as fh:
                    fh.write(phylip([(ref, rs), (r, qs)]))
        if write_multi and cc_cols:
            multi = [(ref, "".join(refc[i] for i in cc_cols))]
            for k, r in enumerate(rows):
                multi.append((r, "".join(fates[k][i][0] for i in cc_cols)))
            with open(os.path.join(out_dir, f"{name}.{tid}.all.phy"), "w") as fh:
                fh.write(phylip(multi))
        with open(os.path.join(out_dir, f"{name}.{tid}.pairs.tsv"), "w") as fh:
            cols = ["informant", "dist", "codons", "retained_pairwise", "retained_complete_case",
                    "written", "identical", "identity", *REASONS]
            fh.write("\t".join(cols) + "\n")
            for e in table:
                fh.write("\t".join("" if e[c] is None else (f"{e[c]:.4f}" if isinstance(e[c], float) else str(e[c]))
                                   for c in cols) + "\n")
        side = {"tool": "pairwise_codons.py", "version": TOOL_VERSION, "stem": name, "reference": ref,
                "chrom": chrom, "window": [start, end], "transcript": tid, "gene": t.get("gene"),
                "strand": strand, "codons": n, "codons_in_cds": codons_in_cds, "cds_length_mod3": mod3,
                "reference_stop_removed": stop_removed, "reference_internal_stops": internal_stops,
                "policy_written": "complete-case" if complete_case else "pairwise",
                "keep_stops": keep_stops, "isoforms": transcript or isoforms,
                "dropped_species": sorted(drop), "informants": rows, "block_selection": block_stats,
                "complete_case_columns": len(cc_cols), "inputs_sha256": checks, "pairs": table}
        with open(os.path.join(out_dir, f"{name}.{tid}.sidecar.json"), "w") as fh:
            json.dump(side, fh, indent=1, sort_keys=True)
        results.append(side)
        if not quiet:
            kept_any = sum(1 for e in table if e["retained_pairwise"])
            print(f"{name}/{tid} ({t.get('gene')}, {strand}): {n} codons, {len(rows)} informants, "
                  f"{kept_any} with >=1 retained codon, complete-case columns {len(cc_cols)}"
                  + (f", cds_length_mod3={mod3}" if mod3 else ""))
    return results


# --------------------------------------------------------------- self-test

def _write_synthetic(d: str) -> str:
    """The three-row example of note 20260909T164106Z-stalin-0016 plus a
    minus-strand spliced gene, an insertion inside a codon, an informant
    stop, and a tree leaf absent from every block."""
    # reference window chrX:1000-1050 (50 bases)
    #  plus gene A: CDS 1000-1015 = ATG GCC TTT CAA TAA (terminal stop removed -> 4 codons)
    #  minus gene B: exons 1030-1036 and 1039-1050; CDS 1030-1036 + 1039-1048 (15 bases)
    #    transcription order = revcomp(ref[1039:1048]) + revcomp(ref[1030:1036])
    #    = GAAGCTATG + TGAAGT -> GAA GCT ATG TGA AGT (internal TGA at index 3)
    ref = "ATGGCCTTTCAATAA" + "G" * 15 + "ACTTCA" + "GGG" + "CATAGCTTCTT"
    assert len(ref) == 50, len(ref)
    tail = ref[12:]
    # one block; a gap column in the reference after position 1004 carries the insertion
    rows = {"ref": "ATGGC-CTTTCAA" + tail,
            "inf": "ATGGC-TTTCCAG" + tail,   # clean: codons ATG GCT TTC CAG
            "thd": "ATGGC-NTT-CAA" + tail,   # N in codon 2, gap in codon 3
            "ins": "ATGGCACTGACAA" + tail}   # insertion inside codon 2, TGA at codon 3
    maf = ["##maf version=1", "a score=0"]
    for i, (nm, sq) in enumerate(rows.items()):
        chrom = "chrX" if nm == "ref" else "chr1"
        maf.append(f"s {nm}.{chrom} {1000 if nm == 'ref' else 5000 * (i + 1)} {len(sq.replace('-', ''))} + 100000 {sq}")
    maf.append("")
    with open(os.path.join(d, "syn.maf"), "w") as fh:
        fh.write("\n".join(maf))
    with open(os.path.join(d, "syn.fa"), "w") as fh:
        fh.write(">chrX:1000-1050\n" + ref + "\n")
    with open(os.path.join(d, "syn.nh"), "w") as fh:
        fh.write("(((ref:0.1,inf:0.2):0.3,(thd:0.4,ins:0.5):0.1):0.5,far:1.0);\n")
    ann = {"transcripts": [
        {"id": "A.1", "gene": "A", "strand": "+", "exons": [[1000, 1015]], "cds": [[1000, 1015]]},
        {"id": "B.1", "gene": "B", "strand": "-", "exons": [[1030, 1036], [1039, 1050]], "cds": [[1030, 1036], [1039, 1048]]},
    ]}
    with open(os.path.join(d, "syn.annotation.json"), "w") as fh:
        json.dump(ann, fh)
    with open(os.path.join(d, "syn.manifest.json"), "w") as fh:
        json.dump({"assembly": "ref", "track": "multiz5way", "window": {"chrom": "chrX", "start": 1000, "end": 1050}}, fh)
    return os.path.join(d, "syn")


def self_test() -> int:
    fails = 0

    def check(cond: bool, what: str) -> None:
        nonlocal fails
        print(("ok   " if cond else "FAIL ") + what)
        if not cond:
            fails += 1

    with tempfile.TemporaryDirectory() as d:
        stem = _write_synthetic(d)
        out = os.path.join(d, "out")
        res = export(stem, out, None, "union", None, set(), False, False, True, "all", quiet=True)
        by = {r["transcript"]: r for r in res}
        a = by["A.1"]
        pa = {e["informant"]: e for e in a["pairs"]}
        check(a["informants"] == ["inf", "thd", "ins", "far"], "rows follow the tree's leaf order, reference excluded")
        check(a["codons"] == 4 and a["reference_stop_removed"], "gene A: 4 codons after removing the terminal TAA")
        check(pa["inf"]["retained_pairwise"] == 4, "clean informant keeps all 4 codons pairwise")
        check(pa["thd"]["retained_pairwise"] == 2 and pa["thd"]["ambiguous"] == 1 and pa["thd"]["gap"] == 1,
              "N and gap informant keeps 2, drops one ambiguous and one gap")
        check(pa["ins"]["insertion"] == 1 and pa["ins"]["informant_stop"] == 1 and pa["ins"]["retained_pairwise"] == 2,
              "insertion inside codon 2 and TGA at codon 3 are dropped for the third informant")
        check(pa["far"]["unaligned"] == 4 and pa["far"]["retained_pairwise"] == 0,
              "a tree leaf in no block is unaligned for every codon")
        check(a["complete_case_columns"] == 0 and not os.path.exists(os.path.join(out, "syn.A.1.far.phy")),
              "complete-case over all rows keeps nothing while an unaligned row is present; no empty PHYLIP")
        with open(os.path.join(out, "syn.A.1.inf.phy")) as fh:
            lines = fh.read().splitlines()
        check(lines[0] == "2 12" and lines[1].split()[1] == "ATGGCCTTTCAA" and lines[2].split()[1] == "ATGGCTTTCCAG",
              "pairwise PHYLIP for the clean informant carries the full CDS pair")
        check(pa["inf"]["identical"] == 1 and abs(pa["inf"]["identity"] - 0.25) < 1e-9,
              "identity counts equal codons (ATG only)")
        check(abs(pa["inf"]["dist"] - 0.3) < 1e-9 and abs(pa["ins"]["dist"] - 1.0) < 1e-9,
              "patristic distances from the tree")
        res1 = export(stem, out, "A.1", "union", ["inf", "thd", "ins"], set(), False, False, True, "all", quiet=True)
        check(res1[0]["complete_case_columns"] == 2 and all(e["retained_complete_case"] == 2 for e in res1[0]["pairs"]),
              "complete-case over the three aligned rows keeps codons 1 and 4 (stalin's [1, 4])")
        with open(os.path.join(out, "syn.A.1.all.phy")) as fh:
            m = fh.read().splitlines()
        check(m[0] == "4 6" and [ln.split()[1] for ln in m[1:]] == ["ATGCAA", "ATGCAG", "ATGCAA", "ATGCAA"],
              "multi PHYLIP holds the two complete-case codons for all four rows")
        b = by["B.1"]
        pb = {e["informant"]: e for e in b["pairs"]}
        check(b["codons_in_cds"] == 5 and b["codons"] == 4 and b["strand"] == "-"
              and b["reference_internal_stops"] == [3] and not b["reference_stop_removed"],
              "minus-strand spliced gene: 5 codons in transcription order, internal TGA excluded and reported")
        with open(os.path.join(out, "syn.B.1.inf.phy")) as fh:
            lb = fh.read().splitlines()
        check(lb[1].split()[1] == "GAAGCTATGAGT" and lb[2].split()[1] == "GAAGCTATGAGT" and pb["inf"]["retained_pairwise"] == 4,
              "minus-strand codons are reverse complemented and span the splice junction")
        check(pb["far"]["unaligned"] == 4 and pb["ins"]["retained_pairwise"] == 4,
              "gene B: absent leaf unaligned for all 4 codons, aligned rows keep all 4")
        # complete-case written pairs and --informants subset
        res2 = export(stem, out, "A.1", "union", ["inf", "thd"], set(), True, False, False, "all", quiet=True)
        p2 = {e["informant"]: e for e in res2[0]["pairs"]}
        check(res2[0]["complete_case_columns"] == 2 and p2["inf"]["written"] == 2 and p2["inf"]["retained_pairwise"] == 4,
              "--complete-case over {inf, thd} writes stalin's [1,4]: 2 of 4 codons for the clean pair")
        with open(os.path.join(out, "syn.A.1.inf.phy")) as fh:
            l2 = fh.read().splitlines()
        check(l2[0] == "2 6" and l2[1].split()[1] == "ATGCAA", "written pair is the complete-case columns")
        res3 = export(stem, out, "A.1", "union", ["ins"], set(), False, True, False, "all", quiet=True)
        check(res3[0]["pairs"][0]["retained_pairwise"] == 3 and res3[0]["pairs"][0]["informant_stop"] == 0,
              "--keep-stops retains the informant TGA codon")
    print(f"{'PASS' if not fails else 'FAIL'}: {17 - fails}/17 checks")
    return 1 if fails else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--stem", help="path prefix of a fetch_window.py output set")
    ap.add_argument("--out", help="output directory")
    ap.add_argument("--transcript", help="transcript id to export (default: --isoforms policy over all loci)")
    ap.add_argument("--isoforms", default="longest-cds", choices=["union", "longest-cds"],
                    help="which isoforms to export when --transcript is not given (default longest-cds)")
    ap.add_argument("--informants", help="comma-separated informant names to keep (default all rows)")
    ap.add_argument("--drop-species", default="", help="comma-separated rows to drop (leakage rule)")
    ap.add_argument("--complete-case", action="store_true",
                    help="write the columns retained across all rows instead of per pair")
    ap.add_argument("--keep-stops", action="store_true", help="keep informant codons that are stops")
    ap.add_argument("--write-multi", action="store_true", help="also write the complete-case multi-row PHYLIP")
    ap.add_argument("--transcript-types", default="benchmark")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if not a.stem or not a.out:
        ap.error("--stem and --out are required")
    drop = {s.strip() for s in a.drop_species.split(",") if s.strip()}
    informants = [s.strip() for s in a.informants.split(",") if s.strip()] if a.informants else None
    export(a.stem, a.out, a.transcript, a.isoforms, informants, drop, a.complete_case, a.keep_stops,
           a.write_multi, a.transcript_types, a.quiet)
    return 0


if __name__ == "__main__":
    sys.exit(main())
