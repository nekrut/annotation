"""Expanded reference recurrence for candidate A's grammar (proposal 3.1, 3.2).

This is the slow specification oracle: every intron carries its m-1
mandatory positions as explicit states, and every state is a Python tuple
in a dictionary per boundary. Nothing here is meant to be fast; the
delayed-entry decoder is checked against it on small lattices.

Coordinates are oriented and zero-based. A state at boundary t describes
bases strictly before t; a transition consumes x[t]. States:

    ('U',)                 outside a coding chain
    ('S', q)               initiator not yet complete, q of length 1 or 2
    ('E', q)               ordinary CDS, q of length 0, 1 or 2 (spliced prefix)
    ('I', c, k)            k-th mandatory intronic base consumed, 1 <= k < m, c in E or S
    ('T', c, r)            intron tail component r, at least m bases consumed

Both forward (log-sum-exp) and Viterbi (max) run the same transition
generator. Viterbi picks the best joint chain and component assignment, as
section 3.2 requires. Ambiguous observed bases branch over the permitted
bases with a uniform log prior, summed in forward and maximized in Viterbi.
Edge partials (3.1's explicit partial-end prior) are not implemented yet:
paths start and end in U.
"""
from dataclasses import dataclass, field
from math import log, exp, inf, isinf
from typing import Dict, List, Optional, Tuple

from .codes import GeneticCode, TABLES, permitted_bases
from .scores import Scores

U = ("U",)
NEG = -inf


def logsumexp(vals):
    m = max(vals)
    if isinf(m):
        return m
    return m + log(sum(exp(v - m) for v in vals))


@dataclass(frozen=True)
class DurationMixture:
    """Pr(l | p) = sum_r pi[p][r] (1-q[p][r]) q[p][r]**(l-m), l >= m (proposal 3.2)."""
    m: int = 20
    pi: Tuple[Tuple[float, ...], ...] = ((1.0,), (1.0,), (1.0,))     # pi[p][r]
    q: Tuple[Tuple[float, ...], ...] = ((0.5,), (0.5,), (0.5,))       # q[p][r]

    def __post_init__(self):
        if self.m < 1:
            raise ValueError("m must be at least 1")
        for p in range(3):
            if len(self.pi[p]) != len(self.q[p]):
                raise ValueError("pi and q need the same number of components")
            if abs(sum(self.pi[p]) - 1.0) > 1e-9:
                raise ValueError("pi[p] must sum to one")
            if not all(0.0 < x < 1.0 for x in self.q[p]):
                raise ValueError("0 < q < 1 required")

    @property
    def R(self):
        return len(self.pi[0])

    def log_pi(self, p, r):
        return log(self.pi[p][r])

    def log_q(self, p, r):
        return log(self.q[p][r])

    def log_1mq(self, p, r):
        return log(1.0 - self.q[p][r])

    def log_prob(self, p, length):
        """Closed-form log Pr(length | p); -inf below m."""
        if length < self.m:
            return NEG
        return logsumexp([self.log_pi(p, r) + self.log_1mq(p, r) + (length - self.m) * self.log_q(p, r)
                          for r in range(self.R)])


@dataclass
class Segment:
    kind: str            # 'cds' or 'intron'
    start: int           # oriented, zero-based, inclusive
    end: int             # exclusive
    prefix_len: int = 0  # CDS bases emitted mod 3 before this segment (p)
    component: Optional[int] = None   # duration component for introns

    @property
    def gff_phase(self):
        """GFF3 phase of a CDS row beginning after p emitted bases: (3-p) % 3."""
        return (3 - self.prefix_len) % 3


@dataclass
class Chain:
    segments: List[Segment] = field(default_factory=list)
    partial_5: bool = False
    partial_3: bool = False
    uncertain: bool = False   # an ambiguous base was consumed as CDS

    def cds(self):
        return [s for s in self.segments if s.kind == "cds"]

    def introns(self):
        return [s for s in self.segments if s.kind == "intron"]

    def spliced(self, x):
        return "".join(x[s.start:s.end] for s in self.cds())


class ReferenceDecoder:
    def __init__(self, code: GeneticCode = TABLES[1], duration: DurationMixture = DurationMixture()):
        self.code = code
        self.dur = duration
        self._prefixes = code.initiator_prefixes()

    # -- coding transition of state c consuming one concrete base b at boundary t
    def _coding_step(self, c, b, t, sc: Scores):
        """Yield (next_state, score) for the ordinary coding transition of c on
        base b. The next state may be U (chain terminated by a stop)."""
        kind, q = c
        p = len(q)
        codon = q + b
        emit = sc.cds[p][t]
        if kind == "S":
            if len(codon) < 3:
                if codon in self._prefixes:
                    yield ("S", codon), emit
            elif self.code.is_initiator(codon):
                yield ("E", ""), emit
        elif kind == "E":
            if len(codon) < 3:
                yield ("E", codon), emit
            elif self.code.is_stop(codon):
                yield U, emit + sc.stop[t]
            else:
                yield ("E", ""), emit

    def _coding_transitions(self, c, x, t, sc):
        """Coding transitions of c consuming observed symbol x[t], branching
        over permitted bases with a uniform prior; yields (next, score, base)."""
        bases = permitted_bases(x[t])
        prior = -log(len(bases))
        for b in bases:
            for nxt, s in self._coding_step(c, b, t, sc):
                yield nxt, s + prior, b

    def transitions(self, state, x, t, sc: Scores):
        """All transitions out of `state` consuming x[t]: (next, score, tag).
        tag records the latent base (coding moves) or None."""
        m, R = self.dur.m, self.dur.R
        kind = state[0]
        if kind == "U":
            yield U, sc.u[t], None
            for b in permitted_bases(x[t]):
                if b in self._prefixes:
                    prior = -log(len(permitted_bases(x[t])))
                    yield ("S", b), sc.start[t] + sc.cds[0][t] + prior, b
            return
        if kind in ("S", "E"):
            for nxt, s, b in self._coding_transitions(state, x, t, sc):
                yield nxt, s, b
            p = len(state[1])
            # donor: open an intron without changing the coding state
            if m == 1:
                for r in range(R):
                    yield ("T", state, r), sc.donor[t] + sc.intron[p][t] + self.dur.log_pi(p, r), None
            else:
                yield ("I", state, 1), sc.donor[t] + sc.intron[p][t], None
            return
        if kind == "I":
            _, c, k = state
            p = len(c[1])
            if k + 1 < m:
                yield ("I", c, k + 1), sc.intron[p][t], None
            else:
                for r in range(R):
                    yield ("T", c, r), sc.intron[p][t] + self.dur.log_pi(p, r), None
            return
        if kind == "T":
            _, c, r = state
            p = len(c[1])
            yield state, sc.intron[p][t] + self.dur.log_q(p, r), None
            exit_score = self.dur.log_1mq(p, r) + sc.acceptor[t]
            for nxt, s, b in self._coding_transitions(c, x, t, sc):
                yield nxt, s + exit_score, b
            return
        raise ValueError(f"unknown state {state}")

    def _run(self, x, sc: Scores, viterbi: bool):
        n = len(x)
        if sc.n != n:
            raise ValueError("scores and sequence length differ")
        layers: List[Dict[tuple, float]] = [{U: 0.0}]
        back: List[Dict[tuple, Tuple[tuple, Optional[str]]]] = [{}]
        for t in range(n):
            cur: Dict[tuple, List[float]] = {}
            bp: Dict[tuple, Tuple[float, tuple, Optional[str]]] = {}
            for state, score in layers[t].items():
                if isinf(score):
                    continue
                for nxt, s, tag in self.transitions(state, x, t, sc):
                    tot = score + s
                    if isinf(tot):
                        continue
                    cur.setdefault(nxt, []).append(tot)
                    if viterbi and (nxt not in bp or tot > bp[nxt][0]):
                        bp[nxt] = (tot, state, tag)
            layer = {k: (max(v) if viterbi else logsumexp(v)) for k, v in cur.items()}
            layers.append(layer)
            back.append({k: (v[1], v[2]) for k, v in bp.items()})
        return layers, back

    def partition(self, x, sc: Scores):
        """log of the sum over all legal paths from U at 0 to U at n."""
        layers, _ = self._run(x, sc, viterbi=False)
        return layers[-1].get(U, NEG)

    def viterbi(self, x, sc: Scores):
        """(best score, Chain list) over paths from U at 0 to U at n."""
        layers, back = self._run(x, sc, viterbi=True)
        n = len(x)
        best = layers[n].get(U, NEG)
        if isinf(best):
            return best, []
        # trace states at every boundary n..0
        states = [None] * (n + 1)
        tags = [None] * n
        s = U
        for t in range(n, 0, -1):
            states[t] = s
            prev, tag = back[t][s]
            tags[t - 1] = tag
            s = prev
        states[0] = s
        return best, self._chains(states, tags, x)

    @staticmethod
    def _chains(states, tags, x):
        """Turn the state path into chains of CDS and intron segments."""
        chains: List[Chain] = []
        cur: Optional[Chain] = None
        seg: Optional[Segment] = None
        n = len(tags)
        for t in range(n):
            a, b = states[t], states[t + 1]
            a_kind, b_kind = a[0], b[0]
            # classify the consumed base x[t]
            if a_kind == "U" and b_kind == "U":
                continue
            if a_kind == "U":                       # chain starts: x[t] is the first initiator base
                cur = Chain()
                seg = Segment("cds", t, t + 1, 0)
                if tags[t] is not None and len(permitted_bases(x[t])) > 1:
                    cur.uncertain = True
                continue
            if a_kind in ("S", "E"):
                if b_kind in ("S", "E", "U"):       # coding base
                    if seg is None or seg.kind != "cds":
                        seg = Segment("cds", t, t + 1, len(a[1]))
                    else:
                        seg.end = t + 1
                    if len(permitted_bases(x[t])) > 1:
                        cur.uncertain = True
                    if b_kind == "U":               # terminal stop completed
                        cur.segments.append(seg)
                        chains.append(cur)
                        cur, seg = None, None
                else:                               # donor: x[t] is intronic
                    cur.segments.append(seg)
                    seg = Segment("intron", t, t + 1, len(a[1]))
                continue
            # a is I or T: x[t] intronic unless we exit
            if b_kind in ("I", "T"):
                seg.end = t + 1
                if b_kind == "T":
                    seg.component = b[2]
            else:                                   # exit: x[t] is the first following CDS base
                cur.segments.append(seg)
                c = a[1]
                seg = Segment("cds", t, t + 1, len(c[1]))
                if len(permitted_bases(x[t])) > 1:
                    cur.uncertain = True
                if b_kind == "U":
                    cur.segments.append(seg)
                    chains.append(cur)
                    cur, seg = None, None
        return chains
