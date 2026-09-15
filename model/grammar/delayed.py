"""Delayed-entry recurrence for candidate A's grammar (proposal 3.2).

The optimized form of the reference decoder in `reference.py`. Only the
duration tails T(c, r) are active intron states; the m-1 mandatory
intronic positions are not visited. A donor from coding state c at
boundary s is parked as a pending entry `D[s, c] + d[s]`; exactly m
boundaries later it enters every tail r with

    D[s, c] + d[s] + sum(u[p, j] for j in range(s, s + m)) + log pi[p, r]

where the emission sum is the rolling window of the last m per-phase
intron emissions. Masked (-inf) emissions are not added to the rolling
sum; a per-phase count of masked positions inside the window makes the
jump emission -inf while any is present and lets the finite sum recover
once the last one leaves (engels-0037: a -inf that enters a plain running
sum leaves it as NaN, -inf - -inf, and poisons every later donor). A tail continues with `u[p, t-1] + log q[p, r]` and
exits with `log(1 - q[p, r])`, the acceptor at t and the ordinary coding
transition consuming x[t]. Pending entries younger than m boundaries live
in a ring buffer (a deque of at most m layers). Edge partials use the same
`EdgePrior` as the reference: E0 and J entry at boundary 0, exit at
boundary n from S, E, T, or from a pending donor whose intron is censored
before reaching m bases (the reference's I(c, k) exit).

Traceback returns the same `Chain` objects as the reference: a tail entry
is expanded back into the I(c, 1..m-1) states it replaced, so the two
decoders' Viterbi paths compare exactly.

Checkpoints and seams (3.3): `partition` and `viterbi` accept interior
`seams`; the recurrence is checkpointed at each seam (`Checkpoint`: active
scores, pending donor ring, emission window with finite sum and mask
count) and the next chunk resumes from that object alone. A seam offers
no partial entry or exit, and an intron spanning it keeps its prefix,
component and score; a traceback jump across a seam resumes from the
donor's true boundary in the earlier chunk.

Standard library only; this is the semantics check, not the tensor
implementation. Per boundary it keeps dictionaries keyed by state tuple.
"""
from collections import deque
from dataclasses import dataclass
from math import isinf, isnan
from typing import Dict, List, Optional, Tuple

from .reference import (ReferenceDecoder, DurationMixture, EdgePrior, Chain,
                        U, NEG, logsumexp, prefix_prior)
from .codes import GeneticCode, TABLES, permitted_bases
from .scores import Scores


class DelayedEntryDecoder(ReferenceDecoder):
    """Same interface as ReferenceDecoder; different recurrence."""

    def __init__(self, code: GeneticCode = TABLES[1], duration: DurationMixture = DurationMixture(),
                 edges: Optional[EdgePrior] = None):
        super().__init__(code, duration, edges)

    # Active-state transitions consuming x[t]: coding moves, tail
    # continuation and tail exits. Donors are *not* yielded here; they go
    # to the pending buffer. Yields (next, score, tag).
    def _active(self, state, x, t, sc):
        kind = state[0]
        if kind == "U":
            yield U, sc.u[t], None
            bases = permitted_bases(x[t])
            for b in bases:
                if b in self._prefixes:
                    yield ("S", b), sc.start[t] + sc.cds[0][t] - _log(len(bases)), b
            return
        if kind == "E0":
            for nxt, s, b in self._coding_transitions(("E", state[1]), x, t, sc):
                yield nxt, s, b
            return
        if kind in ("S", "E"):
            for nxt, s, b in self._coding_transitions(state, x, t, sc):
                yield nxt, s, b
            return
        if kind in ("T", "J"):
            _, c, r = state
            p = len(c[1])
            yield state, sc.intron[p][t] + self.dur.log_q(p, r), None
            if kind == "J" and t == 0:
                return                      # residual intron: at least one observed base before closing
            exit_score = self.dur.log_1mq(p, r) + sc.acceptor[t]
            for nxt, s, b in self._coding_transitions(c, x, t, sc):
                yield nxt, s + exit_score, b
            return
        raise ValueError(f"unknown state {state}")

    def _run(self, x, sc: Scores, viterbi: bool, seams=()):
        """Forward (or Viterbi) pass over the whole sequence. `seams` lists
        interior boundaries at which the recurrence is checkpointed and
        resumed from the checkpoint alone (proposal 3.3): the chunk after a
        seam sees nothing of the chunk before it except the `Checkpoint`.
        Interior seams offer no partial entry or exit."""
        n = len(x)
        self._check_input(x, sc)
        seams = sorted(set(seams))
        if any(t <= 0 or t >= n for t in seams):
            raise ValueError("seams must be interior boundaries 0 < t < n")
        layers: List[Dict[tuple, float]] = []
        back: List[Dict[tuple, Tuple[tuple, Optional[str]]]] = []
        state = Checkpoint.initial(self.initial())
        for t1 in seams + [n]:
            state = state.copy()            # the previous chunk's scratch is not reused
            chunk_layers, chunk_back, state = self._run_chunk(x, sc, viterbi, state, t1)
            layers += chunk_layers
            back += chunk_back
        layers.append(state.layer)
        back.append(state.back)             # pointers into boundary n
        # censored donors still pending at boundary n (reference: I(c,k) exit)
        tail_exits = {}
        if self.edges is not None:
            for s, parked in state.pending:
                for c, base in parked.items():
                    p = len(c[1])
                    emit = sum(sc.intron[p][j] for j in range(s, n))
                    tail_exits[("I", c, n - s)] = (base + emit, ("P", s, c))
        return layers, back, tail_exits

    def _run_chunk(self, x, sc: Scores, viterbi: bool, cp: "Checkpoint", t1: int):
        """Advance the recurrence from checkpoint `cp` (at boundary cp.t) to
        boundary t1, consuming x[cp.t:t1]. Returns the layers and back
        pointers for boundaries cp.t .. t1-1 and the checkpoint at t1."""
        m, R = self.dur.m, self.dur.R
        combine = max if viterbi else logsumexp
        layers, back = [], []
        pending, window, rolling, masked = cp.pending, cp.window, cp.rolling, cp.masked
        layer = cp.layer
        for t in range(cp.t, t1):
            layers.append(layer)
            back.append(cp.back if t == cp.t else bp_prev)
            cur: Dict[tuple, List[float]] = {}
            bp: Dict[tuple, Tuple[float, tuple, Optional[str]]] = {}

            def add(nxt, tot, prev, tag):
                if isinf(tot) or isnan(tot):
                    return
                cur.setdefault(nxt, []).append(tot)
                if viterbi and (nxt not in bp or tot > bp[nxt][0]):
                    bp[nxt] = (tot, prev, tag)

            # 1. park donors from every S/E state at boundary t
            parked = {}
            for state, score in layer.items():
                if state[0] in ("S", "E") and not isinf(score + sc.donor[t]):
                    parked[state] = score + sc.donor[t]
            pending.append((t, parked))
            # 2. active transitions consuming x[t]
            for state, score in layer.items():
                if isinf(score):
                    continue
                for nxt, s, tag in self._active(state, x, t, sc):
                    add(nxt, score + s, state, tag)
            # 3. rolling window of intron emissions: after consuming x[t] the
            #    window covers boundaries t-m+1 .. t
            window.append([sc.intron[p][t] for p in range(3)])
            for p in range(3):
                if isinf(sc.intron[p][t]):
                    masked[p] += 1
                else:
                    rolling[p] += sc.intron[p][t]
            if len(window) > m:
                old = window.popleft()
                for p in range(3):
                    if isinf(old[p]):
                        masked[p] -= 1
                    else:
                        rolling[p] -= old[p]
            # 4. the donor parked at s = t+1-m enters the tails at boundary t+1
            if len(pending) == m:
                s, parked = pending.popleft()
                for c, base in parked.items():
                    p = len(c[1])
                    if masked[p]:
                        continue                # a forbidden intronic base inside the mandatory interval
                    for r in range(R):
                        add(("T", c, r), base + rolling[p] + self.dur.log_pi(p, r), ("P", s, c), None)
            layer = {k: combine(v) for k, v in cur.items()}
            bp_prev = {k: (v[1], v[2]) for k, v in bp.items()}
        return layers, back, Checkpoint(t1, layer, bp_prev if t1 > cp.t else cp.back,
                                        pending, window, rolling, masked)

    def _finals(self, layers, tail_exits):
        finals = {}
        for k, v in layers[-1].items():
            tot = v + self.terminal(k)
            if not (isinf(tot) or isnan(tot)):
                finals[k] = tot
        for k, (v, _) in tail_exits.items():
            tot = v + self.terminal(k)
            if not (isinf(tot) or isnan(tot)):
                finals[k] = tot
        return finals

    def partition(self, x, sc: Scores, seams=()):
        layers, _, tail_exits = self._run(x, sc, viterbi=False, seams=seams)
        finals = self._finals(layers, tail_exits)
        return logsumexp(list(finals.values())) if finals else NEG

    def viterbi(self, x, sc: Scores, seams=()):
        layers, back, tail_exits = self._run(x, sc, viterbi=True, seams=seams)
        n = len(x)
        finals = self._finals(layers, tail_exits)
        if not finals:
            return NEG, []
        s = max(finals, key=finals.get)
        best = finals[s]
        states = [None] * (n + 1)
        tags = [None] * n
        t = n
        if s[0] == "I":                     # censored pending donor: expand I states back to the donor
            _, c, k = s
            donor_at = n - k
            for j in range(n, donor_at, -1):
                states[j] = ("I", c, j - donor_at)
            tags[donor_at] = None
            t = donor_at
            s = c
        while t > 0:
            states[t] = s
            prev, tag = back[t][s]
            if prev[0] == "P":              # tail entry: replay the m-1 mandatory I states
                _, donor_at, c = prev
                for j in range(t - 1, donor_at, -1):
                    states[j] = ("I", c, j - donor_at)
                for j in range(donor_at, t):
                    tags[j] = None
                t = donor_at
                s = c
            else:
                tags[t - 1] = tag
                t -= 1
                s = prev
        states[0] = s
        return best, self._chains(states, tags, x)


@dataclass
class Checkpoint:
    """Everything the delayed recurrence needs to resume at boundary t
    (proposal 3.3): the active scores, the back pointers into boundary t,
    the pending donors younger than m boundaries, the last m per-phase
    intron emissions with their finite sum and mask count. Traceback for
    the chunk before t stays with that chunk; a tail entry whose donor lies
    before the seam replays from the donor's own boundary, not from an
    invented seam state."""
    t: int
    layer: Dict[tuple, float]
    back: Dict[tuple, Tuple[tuple, Optional[str]]]
    pending: deque
    window: deque
    rolling: List[float]
    masked: List[int]

    @classmethod
    def initial(cls, layer0):
        return cls(0, layer0, {}, deque(), deque(), [0.0, 0.0, 0.0], [0, 0, 0])

    def copy(self):
        return Checkpoint(self.t, dict(self.layer), dict(self.back),
                          deque((s, dict(parked)) for s, parked in self.pending),
                          deque(list(w) for w in self.window), list(self.rolling), list(self.masked))


def _log(k):
    from math import log
    return log(k)
