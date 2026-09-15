#!/usr/bin/env python3
"""Label admission audit for candidate A's grammar (proposal section 3.6).

Given one species' reference GFF3 (and, optionally, its checksummed FASTA)
this module decides which annotated CDS chains may supervise the ordinary
chain loss, and why the rest are masked. It is the structured loader
contract of the proposal: it does not repair the GFF, does not change the
benchmark denominator, and keeps every reason as an independent flag.

Admission order (3.6, in the order the proposal lists it):

1. Benchmark filters, as `benchmark/score.py` applies them: nuclear
   primary sequences (`select_seqids`) and complete protein-coding
   transcripts (no `pseudo=true`, no non-coding `gene_biotype`).
2. One representative per locus (gene): most CDS bases, then longest
   exonic span, then annotation order (`select_isoforms(..., "longest-cds")`
   in `scripts/data/cut_windows.py`, 4.4).
3. Topology masks (4.4): connected components of overlapping full CDS
   spans on one sequence and strand; every representative in a component
   with more than one locus is masked.
4. Metadata audit on the raw rows: exception and translation-exception
   tags (on CDS rows or on the transcript row), partial-end declarations
   read per declared end in transcriptional orientation, overlapping rows,
   positive gaps shorter than m, invalid or inconsistent phase,
   translation-table conflicts.
5. FASTA audit on the spliced CDS in transcriptional orientation: ambiguous
   bases, coordinates beyond the sequence, initiator, terminal stop,
   in-frame stops and length modulo 3, with the missing prefix of a
   5'-partial chain marginalized rather than trimmed from every exon. The
   FASTA audit runs on every benchmark-accepted transcript, not only the
   representatives, because alternative transcripts contribute auxiliary
   boundary targets and those must be as reliable as the representatives'.

A declared 5'-partial end is admissible only where the chain's own 5' end
touches the sequence edge in transcriptional orientation (position 1 on
the plus strand, the last base on the minus strand), a 3'-partial end
only where its own 3' end does; a chain declaring both must satisfy both.
A range declaration on an interior row is an interior partial: it masks
the chain and makes the junction site at that row's boundary unknown.

Transcript identity is the triple (seqid, strand, transcript id)
throughout; FlyBase and other sources reuse a transcript id on more than
one sequence or strand, and every cache and membership test keys on the
full triple so an unrelated record can never overwrite another's audit.

Each representative gets a status (`admitted` or `masked`) and its reason
flags; masked spans are the full CDS span, including introns. The manifest
is one TSV row per representative, and `summary` carries the counts the
proposal asks for before fitting: all / selected / topology-masked /
metadata-masked / sequence-masked / admitted, each reason and the union,
masked oriented bases, and the auxiliary boundary catalog retained or made
unknown.

Standard library only. Usage:

    python3 -m model.labels.admission --species Saccharomyces_cerevisiae \\
        --gff panel/Saccharomyces_cerevisiae/GCF_000146045.2_R64_genomic.gff.gz \\
        --fasta panel/Saccharomyces_cerevisiae/GCF_000146045.2_R64_genomic.fna.gz \\
        --out-dir model/labels/manifests
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import io
import json
import os
import resource
import sys
import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

TX_TYPES = {"mRNA", "transcript", "V_gene_segment", "C_gene_segment", "D_gene_segment",
            "J_gene_segment", "primary_transcript", "lnc_RNA", "ncRNA", "rRNA", "tRNA",
            "snRNA", "snoRNA", "miRNA", "guide_RNA", "antisense_RNA", "RNase_P_RNA",
            "RNase_MRP_RNA", "telomerase_RNA", "SRP_RNA", "scRNA", "vault_RNA", "Y_RNA"}
STOPS = {1: {"TAA", "TAG", "TGA"}, 6: {"TGA"}}
INITIATORS = {"ATG"}
ALT_INITIATORS = {1: {"TTG", "CTG"}, 6: set()}
COMP = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")

# metadata reasons, then sequence reasons; the manifest lists them in this order
META_REASONS = ["exception", "transl_except", "partial_5", "partial_3", "partial_unlocated",
                "overlapping_rows", "short_gap", "phase_invalid", "phase_inconsistent",
                "table_conflict", "frame_length"]
INFORMATIONAL = {"partial_5", "partial_3"}
SEQ_REASONS = ["coords_out_of_range", "ambiguous_base", "no_initiator", "no_stop",
               "internal_stop", "seq_frame_length"]


def _score_module():
    """benchmark/score.py has no package; load it by path for its filters."""
    path = os.path.join(ROOT, "benchmark", "score.py")
    spec = importlib.util.spec_from_file_location("benchmark_score", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _open(path):
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _attrs(field9: str) -> Dict[str, str]:
    out = {}
    for part in field9.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def md5_file(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@dataclass
class Row:
    row_no: int
    seqid: str
    strand: str
    start: int          # 1-based inclusive, as in the file
    end: int
    phase: str          # raw column 8
    attrs: Dict[str, str]


@dataclass
class Transcript:
    tid: str
    gid: str
    seqid: str
    strand: str
    order: int                      # first appearance in file order
    ttype: str = "CDS-only"
    rows: List[Row] = field(default_factory=list)
    exons: List[Tuple[int, int]] = field(default_factory=list)
    tx_attrs: Dict[str, str] = field(default_factory=dict)

    @property
    def key(self):
        """Full identity: a transcript id may recur on another sequence or strand."""
        return (self.seqid, self.strand, self.tid)

    @property
    def cds_len(self):
        return sum(r.end - r.start + 1 for r in self.rows)

    @property
    def exon_len(self):
        return sum(b - a + 1 for a, b in self.exons) if self.exons else self.cds_len

    @property
    def span(self):
        return min(r.start for r in self.rows), max(r.end for r in self.rows)

    def oriented_rows(self):
        rows = sorted(self.rows, key=lambda r: r.start)
        return rows if self.strand == "+" else rows[::-1]


class _SeqShim:
    """What score.select_seqids needs: seq_len and region_attrs."""
    def __init__(self):
        self.seq_len = {}
        self.region_attrs = {}


def load_gff_rows(path: str):
    """Parse a GFF3 keeping every CDS row with its phase and attributes.
    Returns (transcripts, gene_biotype, seqshim)."""
    shim = _SeqShim()
    gene_biotype: Dict[str, Optional[str]] = {}
    tx: Dict[str, Transcript] = {}
    order = 0
    with _open(path) as fh:
        for row_no, line in enumerate(fh, 1):
            if not line or line[0] == "#":
                if line.startswith("##sequence-region"):
                    f = line.split()
                    if len(f) == 4:
                        shim.seq_len[f[1]] = int(f[3]) - int(f[2]) + 1
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            seqid, ftype, strand = f[0], f[2], f[6]
            if ftype == "region" and f[3] == "1":
                shim.region_attrs[seqid] = f[8]
                shim.seq_len.setdefault(seqid, int(f[4]))
                continue
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                continue
            if ftype in ("gene", "pseudogene"):
                a = _attrs(f[8])
                gid = a.get("ID") or a.get("gene_id")
                if gid:
                    gene_biotype[gid] = (a.get("gene_biotype") or a.get("gene_type") or a.get("biotype")
                                         or ("pseudogene" if ftype == "pseudogene" else None))
            elif ftype == "CDS":
                a = _attrs(f[8])
                tid = (a.get("Parent") or a.get("transcript_id") or a.get("ID") or "%s:%d:%s" % (seqid, start, strand)).split(",")[0]
                key = "%s\x00%s\x00%s" % (seqid, strand, tid)
                t = tx.get(key)
                if t is None:
                    order += 1
                    t = tx[key] = Transcript(tid, tid, seqid, strand, order)
                t.rows.append(Row(row_no, seqid, strand, start, end, f[7], a))
            elif ftype == "exon":
                a = _attrs(f[8])
                tid = (a.get("Parent") or a.get("transcript_id") or "").split(",")[0]
                if tid:
                    key = "%s\x00%s\x00%s" % (seqid, strand, tid)
                    t = tx.get(key)
                    if t is None:
                        order += 1
                        t = tx[key] = Transcript(tid, tid, seqid, strand, order)
                    t.exons.append((start, end))
            elif ftype not in ("stop_codon", "start_codon"):
                a = _attrs(f[8])
                tid = a.get("ID") or a.get("transcript_id")
                if tid:
                    key = "%s\x00%s\x00%s" % (seqid, strand, tid)
                    t = tx.get(key)
                    if t is None:
                        order += 1
                        t = tx[key] = Transcript(tid, tid, seqid, strand, order)
                    t.ttype = ftype
                    t.tx_attrs = a
                    gid = a.get("Parent") or a.get("gene_id")
                    if gid:
                        t.gid = gid.split(",")[0]
    for k, t in list(tx.items()):
        if not t.rows:
            del tx[k]
    return list(tx.values()), gene_biotype, shim


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


@dataclass
class Meta:
    """Result of the metadata audit of one transcript."""
    flags: set
    prefix: int                 # missing codon bases before a 5'-partial chain
    five: bool                  # chain's own 5' end declared partial (or unlocated)
    three: bool                 # chain's own 3' end declared partial (or unlocated)
    unknown_donors: set = field(default_factory=set)      # genomic positions made unknown
    unknown_acceptors: set = field(default_factory=set)   # by interior partial declarations

    def __iter__(self):
        # backwards-compatible unpacking: (flags, prefix, five, three)
        return iter((self.flags, self.prefix, self.five, self.three))

    def __getitem__(self, i):
        return (self.flags, self.prefix, self.five, self.three)[i]


def metadata_flags(t: Transcript, m: int, table: int) -> Meta:
    """Section 3.6 metadata audit of one transcript's raw rows."""
    flags = set()
    rows = t.oriented_rows()
    # exception tags on CDS rows or on the transcript row
    for r in rows:
        if "exception" in r.attrs:
            flags.add("exception")
        if "transl_except" in r.attrs:
            flags.add("transl_except")
    if "exception" in t.tx_attrs:
        flags.add("exception")
    if "transl_except" in t.tx_attrs:
        flags.add("transl_except")
    # partial declarations, read per row in transcriptional orientation: a
    # start_range is the row's 5' boundary on the plus strand and its 3'
    # boundary on the minus strand, and the other way round for end_range
    five = three = False
    unknown_donors, unknown_acceptors = set(), set()
    last = len(rows) - 1
    for i, r in enumerate(rows):
        has_sr, has_er = "start_range" in r.attrs, "end_range" in r.attrs
        five_side, three_side = (has_sr, has_er) if t.strand == "+" else (has_er, has_sr)
        if five_side:
            if i == 0:
                five = True
            else:
                # the acceptor of the intron before this row is uncertain
                flags.add("partial_unlocated")
                unknown_acceptors.add(r.start - 1 if t.strand == "+" else r.end + 1)
        if three_side:
            if i == last:
                three = True
            else:
                # the donor of the intron after this row is uncertain
                flags.add("partial_unlocated")
                unknown_donors.add(r.end + 1 if t.strand == "+" else r.start - 1)
    if not five and not three and not flags & {"partial_unlocated"} and \
            any(r.attrs.get("partial") == "true" for r in t.rows):
        flags.add("partial_unlocated")
        five = three = True
    if five:
        flags.add("partial_5")
    if three:
        flags.add("partial_3")
    # overlapping rows and short gaps
    srt = sorted(t.rows, key=lambda r: r.start)
    for a, b in zip(srt, srt[1:]):
        gap = b.start - a.end - 1
        if gap < 0:
            flags.add("overlapping_rows")
        elif 0 < gap < m:
            flags.add("short_gap")
    # table
    declared = {r.attrs.get("transl_table") for r in t.rows if r.attrs.get("transl_table")}
    if any(d != str(table) for d in declared):
        flags.add("table_conflict")
    # phase
    phases = []
    for r in rows:
        if r.phase not in ("0", "1", "2"):
            flags.add("phase_invalid")
            phases.append(None)
        else:
            phases.append(int(r.phase))
    p = 0
    if "phase_invalid" not in flags:
        h = phases[0]
        if five:
            p = (-h) % 3
        elif h != 0:
            flags.add("phase_inconsistent")
        for r, h in zip(rows, phases):
            if h != (-p) % 3:
                flags.add("phase_inconsistent")
            p = (p + r.end - r.start + 1) % 3
        if not three and p != 0:
            flags.add("frame_length")
    prefix = (-phases[0]) % 3 if five and phases[0] is not None else 0
    return Meta(flags, prefix, five, three, unknown_donors, unknown_acceptors)


def sequence_flags(t: Transcript, seq: str, prefix: int, five: bool, three: bool, table: int):
    """Section 3.6 FASTA audit of one transcript against its sequence.
    `prefix` is the number of missing codon bases before the observed CDS
    of a 5'-partial chain; they are marginalized, never trimmed."""
    flags = set()
    n = len(seq)
    if any(r.start < 1 or r.end > n for r in t.rows):
        flags.add("coords_out_of_range")
        return flags, None
    spliced = "".join(seq[r.start - 1:r.end] for r in sorted(t.rows, key=lambda r: r.start)).upper()
    if t.strand == "-":
        spliced = revcomp(spliced)
    if any(ch not in "ACGT" for ch in spliced):
        flags.add("ambiguous_base")
    stops = STOPS[table]
    if not five and spliced[:3] not in INITIATORS:
        flags.add("no_initiator")
    # a 5'-partial chain: the first 3-prefix observed bases belong to the
    # codon whose `prefix` bases are missing; that codon is marginalized (every
    # observed suffix has a non-stop completion), the rest is read in frame
    body = spliced[3 - prefix:] if prefix else spliced
    codons = [body[i:i + 3] for i in range(0, len(body) - len(body) % 3, 3)]
    if not three:
        if len(body) % 3:
            flags.add("seq_frame_length")
        if not codons or codons[-1] not in stops:
            flags.add("no_stop")
        internal = codons[:-1]
    else:
        internal = codons
    if any(c in stops for c in internal):
        flags.add("internal_stop")
    return flags, spliced


@dataclass
class AuditResult:
    species: str
    manifest_rows: List[dict]
    summary: dict


def _components(reps: List[Transcript]):
    """Connected components of overlapping full CDS spans per (seqid, strand)."""
    masked = set()
    spans_masked = []
    by = defaultdict(list)
    for t in reps:
        s, e = t.span
        by[(t.seqid, t.strand)].append((s, e + 1, t))     # half-open
    n_comp = 0
    for key, items in by.items():
        items.sort(key=lambda x: (x[0], x[1]))
        cur = []
        cur_end = -1
        for s, e, t in items:
            if cur and s >= cur_end:
                if len({u.gid for u in cur}) > 1:
                    n_comp += 1
                    masked.update(u.key for u in cur)
                    spans_masked.append((key, min(u.span[0] for u in cur), max(u.span[1] for u in cur) + 1))
                cur, cur_end = [], -1
            cur.append(t)
            cur_end = max(cur_end, e)
        if cur and len({u.gid for u in cur}) > 1:
            n_comp += 1
            masked.update(u.key for u in cur)
            spans_masked.append((key, min(u.span[0] for u in cur), max(u.span[1] for u in cur) + 1))
    return masked, n_comp, spans_masked


def _iter_fasta(path):
    name, buf = None, []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def audit_species(species: str, gff: str, fasta: Optional[str], m: int = 20, table: int = 1,
                  log=print) -> AuditResult:
    score = _score_module()
    transcripts, gene_biotype, shim = load_gff_rows(gff)
    seq_reasons = {}
    keep_seq = score.select_seqids(shim, None, reasons=seq_reasons)
    counts = {"cds_transcripts_all": len(transcripts)}
    # 1. benchmark filters
    accepted = []
    filtered = defaultdict(int)
    for t in transcripts:
        if t.seqid not in keep_seq:
            filtered["sequence: " + seq_reasons.get(t.seqid, "not selected")] += 1
            continue
        if any(r.attrs.get("pseudo") == "true" for r in t.rows) or t.tx_attrs.get("pseudo") == "true":
            filtered["pseudogene"] += 1
            continue
        bt = gene_biotype.get(t.gid)
        if bt is not None and (bt.endswith("pseudogene")):
            filtered["pseudogene"] += 1
            continue
        if bt is not None and bt != "protein_coding":
            filtered["biotype=" + bt] += 1
            continue
        accepted.append(t)
    counts["benchmark_filtered"] = dict(sorted(filtered.items()))
    counts["accepted_transcripts"] = len(accepted)
    # 2. representatives: most CDS bases, longest exonic span, annotation order
    by_gene = defaultdict(list)
    for t in accepted:
        by_gene[(t.seqid, t.strand, t.gid)].append(t)
    reps = []
    for members in by_gene.values():
        members.sort(key=lambda t: (-t.cds_len, -t.exon_len, t.order))
        reps.append(members[0])
    reps.sort(key=lambda t: t.order)
    counts["loci"] = len(by_gene)
    counts["representatives"] = len(reps)
    counts["nonrepresentative_omitted"] = len(accepted) - len(reps)
    # 3. topology masks
    topo_masked, n_comp, topo_spans = _components(reps)
    counts["topology_components"] = n_comp
    counts["topology_masked"] = len(topo_masked)
    # 4. metadata audit, on every accepted transcript (alternatives feed the
    #    auxiliary catalog); keyed on the full identity
    meta: Dict[tuple, Meta] = {t.key: metadata_flags(t, m, table) for t in accepted}
    # 5. FASTA audit, streaming one sequence at a time, on every accepted transcript
    seqf = {t.key: (set(), None) for t in accepted}
    fasta_seen = set()
    if fasta:
        need = defaultdict(list)
        for t in accepted:
            need[t.seqid].append(t)
        for name, seq in _iter_fasta(fasta):
            if name not in need:
                continue
            fasta_seen.add(name)
            for t in need[name]:
                mt = meta[t.key]
                seqf[t.key] = sequence_flags(t, seq, mt.prefix, mt.five, mt.three, table)
            need.pop(name)
        missing = sorted(need)
        for name in missing:
            for t in need[name]:
                seqf[t.key] = ({"coords_out_of_range"}, None)
        counts["fasta_sequences_missing"] = len(missing)

    def own_edge(t: Transcript, end: str) -> bool:
        """Does the chain's declared `end` ('5' or '3') touch the sequence edge
        in transcriptional orientation?"""
        s, e = t.span
        n = shim.seq_len.get(t.seqid, -1)
        low = (end == "5") == (t.strand == "+")     # 5' end of a plus chain, 3' end of a minus chain
        return s == 1 if low else e == n

    # assemble
    rows_out = []
    reason_counts = defaultdict(int)
    masked_spans = list(topo_spans)
    n_meta = n_seq = n_adm = 0
    for t in reps:
        mt = meta[t.key]
        mflags, p, five, three = mt.flags, mt.prefix, mt.five, mt.three
        sflags, _ = seqf[t.key]
        reasons = []
        if t.key in topo_masked:
            reasons.append("topology")
        reasons += [r for r in META_REASONS if r in mflags]
        reasons += [r for r in SEQ_REASONS if r in sflags]
        # a declared partial end is masked unless that end really touches the
        # sequence edge; each declared end is checked at its own edge
        if (five or three) and "partial_unlocated" not in mflags:
            if (five and not own_edge(t, "5")) or (three and not own_edge(t, "3")):
                reasons.append("partial_away_from_edge")
        for r in reasons:
            reason_counts[r] += 1
        # a declared partial end at the true sequence edge is admissible in the
        # edge-initialized grammar (3.1, 3.6): partial_5/partial_3 are recorded
        # but do not mask by themselves; away from the edge or unlocated they do
        masking = [r for r in reasons if r not in INFORMATIONAL]
        if any(r not in SEQ_REASONS and r != "topology" for r in masking):
            n_meta += 1
        if sflags:
            n_seq += 1
        status = "admitted" if not masking else "masked"
        if status == "admitted":
            n_adm += 1
        else:
            s, e = t.span
            masked_spans.append(((t.seqid, t.strand), s, e + 1))
        s, e = t.span
        rows_out.append({"transcript": t.tid, "gene": t.gid, "type": t.ttype, "seqid": t.seqid,
                         "strand": t.strand, "cds_start": s, "cds_end": e, "rows": len(t.rows),
                         "cds_len": t.cds_len, "missing_prefix": p, "status": status,
                         "reasons": ",".join(reasons) or "-"})
    # auxiliary catalog over all accepted transcripts: distinct sites. The
    # annotated catalog is every declared site; a site is retained when at
    # least one transcript observes it reliably (its end not declared partial
    # or unlocated, its coordinates in range, and its sequence passing the
    # initiator/stop check for starts and stops), otherwise unknown: suppressed
    # from both positive and negative supervision, never negative
    annotated = {"start": set(), "stop": set(), "donor": set(), "acceptor": set()}
    reliable = {"start": set(), "stop": set(), "donor": set(), "acceptor": set()}
    for t in accepted:
        rows = t.oriented_rows()
        mt = meta[t.key]
        sflags = seqf[t.key][0]
        first, last = rows[0], rows[-1]
        start_pos = first.start if t.strand == "+" else first.end
        stop_pos = last.end if t.strand == "+" else last.start
        key = (t.seqid, t.strand)
        in_range = "coords_out_of_range" not in sflags
        unlocated = "partial_unlocated" in mt.flags and not (mt.unknown_donors or mt.unknown_acceptors)
        annotated["start"].add(key + (start_pos,))
        annotated["stop"].add(key + (stop_pos,))
        if in_range and not mt.five and not unlocated and "no_initiator" not in sflags:
            reliable["start"].add(key + (start_pos,))
        if in_range and not mt.three and not unlocated and "no_stop" not in sflags:
            reliable["stop"].add(key + (stop_pos,))
        srt = sorted(t.rows, key=lambda r: r.start)
        for a, b in zip(srt, srt[1:]):
            gap = b.start - a.end - 1
            if gap < m:
                continue                                      # short gaps are not junctions (4.4)
            d, ac = (a.end + 1, b.start - 1) if t.strand == "+" else (b.start - 1, a.end + 1)
            annotated["donor"].add(key + (d,))
            annotated["acceptor"].add(key + (ac,))
            if in_range and d not in mt.unknown_donors:
                reliable["donor"].add(key + (d,))
            if in_range and ac not in mt.unknown_acceptors:
                reliable["acceptor"].add(key + (ac,))
    # masked oriented bases: union of masked spans per (seqid, strand)
    masked_bases = 0
    by = defaultdict(list)
    for key, s, e in masked_spans:
        by[key].append((s, e))
    for key, ivs in by.items():
        ivs.sort()
        cs, ce = None, None
        for s, e in ivs:
            if cs is None or s > ce:
                if cs is not None:
                    masked_bases += ce - cs
                cs, ce = s, e
            else:
                ce = max(ce, e)
        if cs is not None:
            masked_bases += ce - cs
    counts.update({
        "metadata_masked": n_meta, "sequence_masked": n_seq,
        "masked_union": len(reps) - n_adm, "admitted": n_adm,
        "reasons": dict(sorted(reason_counts.items())),
        "masked_oriented_bases": masked_bases,
        "fasta_checked": bool(fasta),
        "auxiliary_sites_annotated": {k: len(v) for k, v in annotated.items()},
        "auxiliary_sites_retained": {k: len(v) for k, v in reliable.items()},
        "auxiliary_sites_unknown": {k: len(v - reliable[k]) for k, v in annotated.items()},
        "m": m, "table": table,
    })
    return AuditResult(species, rows_out, counts)


COLUMNS = ["transcript", "gene", "type", "seqid", "strand", "cds_start", "cds_end", "rows",
           "cds_len", "missing_prefix", "status", "reasons"]


def write_manifest(res: AuditResult, out_dir: str, gff: str, fasta: Optional[str]):
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, res.species + ".manifest.tsv.gz")
    with gzip.open(path, "wt", encoding="utf-8", compresslevel=9) as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in res.manifest_rows:
            fh.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")
    summary = dict(res.summary)
    summary["species"] = res.species
    summary["gff"] = os.path.basename(gff)
    summary["gff_md5"] = md5_file(gff)
    summary["fasta"] = os.path.basename(fasta) if fasta else None
    summary["fasta_md5"] = md5_file(fasta) if fasta else None
    summary["manifest"] = os.path.basename(path)
    summary["manifest_md5"] = md5_file(path)
    spath = os.path.join(out_dir, res.species + ".summary.json")
    with open(spath, "w") as fh:
        json.dump(summary, fh, indent=1, sort_keys=True)
    return path, spath, summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--m", type=int, default=20)
    ap.add_argument("--table", type=int, default=1)
    a = ap.parse_args(argv)
    t0 = time.time()
    res = audit_species(a.species, a.gff, a.fasta, a.m, a.table)
    res.summary["wall_seconds"] = round(time.time() - t0, 1)
    res.summary["peak_rss_mb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
    path, spath, summary = write_manifest(res, a.out_dir, a.gff, a.fasta)
    print(json.dumps(summary, indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
