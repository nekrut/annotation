#!/usr/bin/env python3
"""Pairwise KA/KS by the Nei and Gojobori (1986) pathway method, with the
Jukes-Cantor correction and the one-sided z-test that asks whether
nonsynonymous substitutions are rarer than synonymous ones.

Standard library only (Python 3.11).  This is the estimator behind the
KA/KS comparative baseline of T-human-010, which reproduces the test of
Nekrutenko, Makova and Li (2002, Genome Research 12:198-202,
doi:10.1101/gr.200901): in a protein-coding region synonymous changes
(KS) outnumber nonsynonymous ones (KA), so a window in which KA/KS is
significantly below one is called coding.  The original study fitted the
ratio by maximum likelihood in codeml and tested it with a likelihood-ratio
statistic (relay note 20260909T164106Z-stalin-0016); codeml is a C
dependency outside the charter's stack, so this module uses the counting
method instead, which is the estimator the same authors' 2003 follow-up
compared it against and the one every textbook z-test of KA < KS uses.
The decision rule (ratio below one *and* a significant difference) is the
same.

Definitions
-----------
For a codon, the synonymous site count at each of its three positions is
the fraction of the three possible single-base changes that leave the amino
acid unchanged; changes to a stop codon count as nonsynonymous, and a stop
codon itself has no sites (it is skipped by the callers).  ``S`` for a pair
is the mean over the two sequences of the summed synonymous sites, ``N`` is
``3L - S``.  Differences between two codons are attributed along the
shortest mutational pathways; when codons differ at two or three positions
every pathway that does not pass through a stop codon is weighted equally
(all pathways if every one passes through a stop).  ``pS = Sd / S`` and
``pN = Nd / N`` are corrected with Jukes-Cantor, ``d = -3/4 ln(1 - 4p/3)``,
with variance ``p (1 - p) / (n (1 - 4p/3)^2)`` (Nei and Kumar 2000, eq.
4.9 and 4.10 in the counting-method section; Kimura 1983 for the variance of
the JC distance).  The z statistic is ``(dS - dN) / sqrt(V(dS) + V(dN))`` and
``p`` is the upper tail of the standard normal (one-sided, H1: dN < dS).

A pair is *uninformative* rather than tested when it has fewer codons than
the caller's floor, when ``pS >= 3/4`` (synonymous saturation; dS is
undefined), or when there is no difference at all (``Sd = Nd = 0``: the
sequences are too close for the test to see selection either way).

Genetic codes: the standard code (translation table 1) by default;
``table=6`` (ciliate: TAA and TAG read as glutamine, TGA stops) is provided
for *Tetrahymena*.  Any other table can be passed as a dict.
"""
from __future__ import annotations

import itertools
import math

BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
TABLE1 = dict(zip(CODONS, AMINO))
TABLE6 = dict(TABLE1, TAA="Q", TAG="Q")
TABLES = {1: TABLE1, 6: TABLE6}


def genetic_code(table) -> dict:
    if isinstance(table, dict):
        return table
    return TABLES[int(table)]


def syn_sites(codon: str, code: dict) -> float:
    """Synonymous sites of one codon (0 to 3); stop codons return 0."""
    aa = code[codon]
    if aa == "*":
        return 0.0
    s = 0.0
    for i in range(3):
        for b in BASES:
            if b == codon[i]:
                continue
            alt = codon[:i] + b + codon[i + 1:]
            if code[alt] == aa:
                s += 1.0 / 3.0
    return s


def pair_differences(c1: str, c2: str, code: dict) -> tuple[float, float]:
    """(synonymous, nonsynonymous) differences between two sense codons,
    averaged over equally weighted shortest pathways that avoid stops."""
    diff = [i for i in range(3) if c1[i] != c2[i]]
    if not diff:
        return 0.0, 0.0
    if len(diff) == 1:
        return (1.0, 0.0) if code[c1] == code[c2] else (0.0, 1.0)
    paths = []
    for order in itertools.permutations(diff):
        cur = c1
        sd = nd = 0
        via_stop = False
        for i in order:
            nxt = cur[:i] + c2[i] + cur[i + 1:]
            if code[nxt] == "*":
                via_stop = True
            if code[cur] == code[nxt]:
                sd += 1
            else:
                nd += 1
            cur = nxt
        paths.append((via_stop, sd, nd))
    clean = [p for p in paths if not p[0]] or paths
    return (sum(p[1] for p in clean) / len(clean), sum(p[2] for p in clean) / len(clean))


class CodonTables:
    """Precomputed site and difference tables for one genetic code."""

    def __init__(self, table=1):
        self.code = genetic_code(table)
        self.stops = {c for c, a in self.code.items() if a == "*"}
        self.sense = [c for c in CODONS if c not in self.stops]
        self.sites = {c: syn_sites(c, self.code) for c in self.sense}
        self.diff = {}
        for a in self.sense:
            for b in self.sense:
                self.diff[(a, b)] = pair_differences(a, b, self.code)


_TABLES: dict = {}


def tables(table=1) -> CodonTables:
    key = table if not isinstance(table, dict) else id(table)
    if key not in _TABLES:
        _TABLES[key] = CodonTables(table)
    return _TABLES[key]


def jukes_cantor(p: float, n: float) -> tuple[float, float]:
    """(distance, variance); (inf, inf) at or beyond saturation."""
    if p >= 0.75:
        return math.inf, math.inf
    if p <= 0.0:
        return 0.0, 0.0
    x = 1.0 - 4.0 * p / 3.0
    d = -0.75 * math.log(x)
    v = p * (1.0 - p) / (n * x * x)
    return d, v


def normal_sf(z: float) -> float:
    """Upper tail of the standard normal."""
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def kaks(pairs, table=1, min_codons: int = 1) -> dict:
    """KA/KS for a list of (reference codon, informant codon) pairs.

    Pairs containing a stop or a non-ACGT letter are dropped and counted.
    Returns a dict with the counts, ``S``, ``N``, ``Sd``, ``Nd``, ``pS``,
    ``pN``, ``dS``, ``dN``, ``ratio`` (dN/dS), ``z``, ``p`` and ``status``:
    ``tested``, ``too_few_codons``, ``identical`` or ``saturated``."""
    T = tables(table)
    L = 0
    dropped = {"stop": 0, "ambiguous": 0}
    S1 = S2 = 0.0
    Sd = Nd = 0.0
    for a, b in pairs:
        a = a.upper()
        b = b.upper()
        if a not in T.code or b not in T.code:
            dropped["ambiguous"] += 1
            continue
        if a in T.stops or b in T.stops:
            dropped["stop"] += 1
            continue
        L += 1
        S1 += T.sites[a]
        S2 += T.sites[b]
        sd, nd = T.diff[(a, b)]
        Sd += sd
        Nd += nd
    S = 0.5 * (S1 + S2)
    N = 3.0 * L - S
    out = {"codons": L, "dropped": dropped, "S": S, "N": N, "Sd": Sd, "Nd": Nd,
           "pS": None, "pN": None, "dS": None, "dN": None, "ratio": None, "z": None, "p": None}
    if L < min_codons or S <= 0 or N <= 0:
        out["status"] = "too_few_codons"
        return out
    pS, pN = Sd / S, Nd / N
    out["pS"], out["pN"] = pS, pN
    if Sd == 0 and Nd == 0:
        out["status"] = "identical"
        out["dS"] = out["dN"] = 0.0
        return out
    dS, vS = jukes_cantor(pS, S)
    dN, vN = jukes_cantor(pN, N)
    out["dS"], out["dN"] = dS, dN
    if math.isinf(dS) or math.isinf(dN):
        out["status"] = "saturated"
        return out
    out["ratio"] = (dN / dS) if dS > 0 else math.inf
    se = math.sqrt(vS + vN)
    if se > 0:
        z = (dS - dN) / se
        out["z"], out["p"] = z, normal_sf(z)
    else:
        out["z"], out["p"] = 0.0, 1.0
    out["status"] = "tested"
    return out


def is_coding(result: dict, alpha: float = 0.05) -> bool:
    """The decision rule: tested, dN/dS < 1, and dN significantly below dS."""
    return result.get("status") == "tested" and result["ratio"] < 1.0 and result["p"] < alpha


# ------------------------------------------------------------------ tests

def self_test() -> int:
    checks = 0
    T = tables(1)
    # 1: site counts of textbook codons (Nei and Kumar 2000, table 4.1 style)
    assert abs(T.sites["TTT"] - 1 / 3) < 1e-9, T.sites["TTT"]        # Phe: third position 1/3
    assert abs(T.sites["CTT"] - 1.0) < 1e-9, T.sites["CTT"]           # Leu CTN: third position fully syn
    assert abs(T.sites["CTG"] - 4 / 3) < 1e-9, T.sites["CTG"]         # Leu CTG: + TTG at position 1
    assert abs(T.sites["ATG"]) < 1e-9 and abs(T.sites["TGG"]) < 1e-9  # Met, Trp
    assert abs(T.sites["AGA"] - 2 / 3) < 1e-9, T.sites["AGA"]         # Arg AGA: CGA (pos 1) and AGG (pos 3)
    checks += 5
    # 2: total synonymous sites over the 61 sense codons of the standard code
    # (61 x 9 = 549 single-base changes: 134 synonymous, 392 nonsynonymous, 23 to a stop)
    total = sum(T.sites.values())
    assert abs(total - 134 / 3) < 1e-9, total
    checks += 1
    # 3: differences along pathways
    assert T.diff[("TTT", "TTC")] == (1.0, 0.0)
    assert T.diff[("TTT", "TTA")] == (0.0, 1.0)
    # TTT -> GTA: via GTT (Val, then GTA syn) = (1 syn, 1 nonsyn); via TTA (Leu nonsyn, then GTA nonsyn) = (0, 2)
    sd, nd = T.diff[("TTT", "GTA")]
    assert abs(sd - 0.5) < 1e-9 and abs(nd - 1.5) < 1e-9, (sd, nd)
    # a pathway through a stop is avoided: TAT (Tyr) -> TGG (Trp): via TGT (Cys) ok, via TAG (stop) excluded
    sd, nd = T.diff[("TAT", "TGG")]
    assert (sd, nd) == (0.0, 2.0), (sd, nd)
    checks += 4
    # 4: Jukes-Cantor
    d, v = jukes_cantor(0.1, 100)
    assert abs(d - 0.107326) < 1e-5 and v > 0, (d, v)
    assert jukes_cantor(0.8, 10)[0] == math.inf
    checks += 2
    # 5: a strongly purifying pair: synonymous third-position changes only
    ref = ["CTT", "GCT", "AAA", "GGT", "TTT", "CCT", "GAA", "AGA", "TCT", "ACT"] * 6
    alt = [c[:2] + ({"T": "C", "A": "G", "C": "T", "G": "A"}[c[2]] if i % 2 == 0 else c[2]) for i, c in enumerate(ref)]
    alt = [c if tables(1).code[c] == tables(1).code[r] else r for c, r in zip(alt, ref)]
    r = kaks(list(zip(ref, alt)))
    assert r["status"] == "tested" and r["Nd"] == 0 and r["Sd"] > 0 and r["ratio"] == 0.0 and r["p"] < 0.01, r
    assert is_coding(r)
    checks += 2
    # 6: a neutral pair: random changes at every position give ratio near 1 and no call
    import random
    rng = random.Random(7)
    ref = [rng.choice(T.sense) for _ in range(300)]
    alt = []
    for c in ref:
        while True:
            m = "".join(rng.choice(BASES) if rng.random() < 0.15 else b for b in c)
            if m not in T.stops:
                break
        alt.append(m)
    r = kaks(list(zip(ref, alt)))
    assert r["status"] == "tested" and 0.6 < r["ratio"] < 1.6, r["ratio"]
    assert not is_coding(r, alpha=0.001), r
    checks += 2
    # 7: identical, too few, stops and ambiguity
    assert kaks([("CTT", "CTT")] * 30)["status"] == "identical"
    assert kaks([("ATG", "ATG")] * 30)["status"] == "too_few_codons"   # no synonymous site at all
    assert kaks([("CTT", "CTT")] * 5, min_codons=10)["status"] == "too_few_codons"
    r = kaks([("ATG", "TAA"), ("CCN", "CCC"), ("CCT", "CCC")] * 10)
    assert r["dropped"] == {"stop": 10, "ambiguous": 10} and r["codons"] == 10, r
    checks += 3
    # 8: saturation
    r = kaks([("TTT", "TTC"), ("TTC", "TTT"), ("CTT", "CTC"), ("CTA", "CTG")] * 10)
    assert r["status"] == "saturated", r["status"]
    checks += 1
    # 9: ciliate code: TAA is a sense codon and has sites
    T6 = tables(6)
    assert "TAA" in T6.sense and "TGA" in T6.stops and T6.code["TAG"] == "Q"
    checks += 1
    # 10: symmetry of the counting in the pair
    a = kaks([("TTT", "GTA")] * 5 + [("CTT", "CTC")] * 5 + [("GCT", "GCT")] * 40)
    b = kaks([("GTA", "TTT")] * 5 + [("CTC", "CTT")] * 5 + [("GCT", "GCT")] * 40)
    assert a["status"] == "tested" and abs(a["dN"] - b["dN"]) < 1e-12 and abs(a["dS"] - b["dS"]) < 1e-12, (a, b)
    checks += 1
    print(f"self-test passed ({checks} checks)")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(self_test())
