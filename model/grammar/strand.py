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
    kind: str            # 'cds' or 'intron'
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


def gff3_rows(chain: Chain, seqid: str, n: int, strand: str, parent: str) -> List[str]:
    """CDS rows (GFF3, one-based inclusive) for one chain; introns are
    implied by the gaps. Partial ends are recorded as attributes so a
    consumer can tell a censored codon from a complete one."""
    rows = []
    for f in genomic_features(chain, n, strand):
        if f.kind != "cds":
            continue
        attrs = [f"Parent={parent}"]
        if chain.uncertain:
            attrs.append("sequence_uncertain=true")
        if chain.partial_5:
            attrs.append("partial_5=true")
        if chain.partial_3:
            attrs.append("partial_3=true")
        rows.append("\t".join([seqid, "grammar", "CDS", str(f.start + 1), str(f.end), ".",
                               strand, str(f.phase), ";".join(attrs)]))
    return rows
