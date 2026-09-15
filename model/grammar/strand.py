"""Strand handling: reverse complement and the oriented-to-genomic map.

The decoders work in oriented, zero-based coordinates (proposal 3.1): the
plus strand is the genomic sequence itself, the minus strand is its reverse
complement, and a chain found on the minus strand of a window of length n
occupies oriented [s, e) exactly where the genomic interval [n-e, n-s)
lies. GFF3 phase is a property of the chain (bases emitted before the
segment, mod 3), so it does not change under the map; only coordinates and
the strand column do. Partial flags keep their chain meaning: `partial_5`
is the chain's 5' end, which on the minus strand is the genomic high end.
"""
from dataclasses import dataclass
from typing import List, Optional

from .reference import Chain

_COMPLEMENT = str.maketrans("ACGTUacgtuRYSWKMBDHVNXryswkmbdhvnx",
                            "TGCAAtgcaaYRSWMKVHDBNXyrswmkvhdbnx")


def reverse_complement(x: str) -> str:
    """Reverse complement with IUPAC ambiguity codes preserved (R<->Y, K<->M,
    B<->V, D<->H; S, W, N, X are self-complementary)."""
    return x.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Feature:
    kind: str            # 'cds', 'intron', 'start_codon' or 'stop_codon'
    start: int           # genomic, zero-based, inclusive
    end: int             # genomic, exclusive
    strand: str          # '+' or '-'
    phase: Optional[int] # GFF3 phase for CDS, None for introns
    component: Optional[int]
    prefix_len: int      # spliced CDS bases mod 3 before this segment along the chain


def genomic_features(chain: Chain, n: int, strand: str) -> List[Feature]:
    """Map a chain decoded on an oriented sequence of length n to genomic
    features on `strand`, sorted by genomic start."""
    if strand not in ("+", "-"):
        raise ValueError("strand must be '+' or '-'")
    out = []
    for s in chain.segments:
        if strand == "+":
            a, b = s.start, s.end
        else:
            a, b = n - s.end, n - s.start
        out.append(Feature(s.kind, a, b, strand,
                           s.gff_phase if s.kind == "cds" else None,
                           s.component, s.prefix_len))
    out.sort(key=lambda f: (f.start, f.end))
    return out


def _codon_pieces(chain: Chain, which: str):
    """Oriented pieces (start, end, bases of the codon before the piece) of
    the chain's initiator (`which` == 'start') or terminal stop ('stop'),
    split across introns when the codon is. A censored end (partial_5 for
    the start, partial_3 for the stop) has no complete codon and yields
    nothing; the terminal stop is part of the CDS and stays there."""
    cds = chain.cds()
    if which == "start":
        if chain.partial_5:
            return []
        segs = cds
    else:
        if chain.partial_3:
            return []
        segs = cds[::-1]
    total = sum(s.end - s.start for s in cds)
    if total < 3:
        return []
    pieces, need = [], 3
    for s in segs:
        take = min(need, s.end - s.start)
        if which == "start":
            pieces.append((s.start, s.start + take, 3 - need))
        else:
            pieces.append((s.end - take, s.end, 0))
        need -= take
        if need == 0:
            break
    if which == "stop":
        # walked from the 3' end: reverse into chain order and set the
        # number of codon bases before each piece
        pieces = pieces[::-1]
        before, out = 0, []
        for a, b, _ in pieces:
            out.append((a, b, before))
            before += b - a
        pieces = out
    return pieces


def codon_features(chain: Chain, n: int, strand: str) -> List[Feature]:
    """Explicit `start_codon` and `stop_codon` features (proposal 3.1) of a
    chain decoded on an oriented sequence of length n, mapped to genomic
    coordinates on `strand`; a codon split across an intron becomes one
    feature per piece with the GFF3 phase of its continuation."""
    if strand not in ("+", "-"):
        raise ValueError("strand must be '+' or '-'")
    out = []
    for which, kind in (("start", "start_codon"), ("stop", "stop_codon")):
        for a, b, before in _codon_pieces(chain, which):
            if strand == "+":
                ga, gb = a, b
            else:
                ga, gb = n - b, n - a
            out.append(Feature(kind, ga, gb, strand, (3 - before) % 3, None, before))
    out.sort(key=lambda f: (f.start, f.end))
    return out


def gff3_rows(chain: Chain, seqid: str, n: int, strand: str, parent: str) -> List[str]:
    """GFF3 rows (one-based inclusive) for one chain: CDS rows, then the
    explicit start_codon and stop_codon features of proposal 3.1, split
    across introns when the codon is. Introns are implied by the gaps; the
    terminal stop stays inside the CDS. Partial ends are recorded as
    attributes so a consumer can tell a censored codon from a complete one,
    and a censored end emits no codon feature."""
    attrs = [f"Parent={parent}"]
    if chain.uncertain:
        attrs.append("sequence_uncertain=true")
    if chain.partial_5:
        attrs.append("partial_5=true")
    if chain.partial_3:
        attrs.append("partial_3=true")
    rows = []
    for f in genomic_features(chain, n, strand):
        if f.kind != "cds":
            continue
        rows.append("\t".join([seqid, "grammar", "CDS", str(f.start + 1), str(f.end), ".",
                               strand, str(f.phase), ";".join(attrs)]))
    for f in codon_features(chain, n, strand):
        rows.append("\t".join([seqid, "grammar", f.kind, str(f.start + 1), str(f.end), ".",
                               strand, str(f.phase), f"Parent={parent}"]))
    return rows
