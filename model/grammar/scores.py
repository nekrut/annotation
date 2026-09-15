"""Per-base score channels the decoder consumes (proposal section 3.1).

Eleven channels per oriented base at boundary t, consuming x[t]:
U (one), CDS by prefix length p (three), intron by p (three), start, stop,
donor, acceptor. Prefix states share the phase emission. A `Scores` of
zeros gives the pure grammar-plus-duration score, which is what the
hand-checkable cases in section 3.4 test.
"""
from dataclasses import dataclass, field
from typing import List


def _zeros(n):
    return [0.0] * n


@dataclass
class Scores:
    n: int
    u: List[float] = field(default=None)
    cds: List[List[float]] = field(default=None)      # cds[p][t]
    intron: List[List[float]] = field(default=None)   # intron[p][t]
    start: List[float] = field(default=None)          # attached to the first initiator base
    stop: List[float] = field(default=None)           # attached to completion of the terminal codon
    donor: List[float] = field(default=None)          # intron opens at boundary t (x[t] is intronic)
    acceptor: List[float] = field(default=None)       # intron closes at boundary t (x[t] is CDS)

    def __post_init__(self):
        n = self.n
        if self.u is None: self.u = _zeros(n)
        if self.cds is None: self.cds = [_zeros(n) for _ in range(3)]
        if self.intron is None: self.intron = [_zeros(n) for _ in range(3)]
        for name in ("start", "stop", "donor", "acceptor"):
            if getattr(self, name) is None:
                setattr(self, name, _zeros(n))
        for arr in [self.u, self.start, self.stop, self.donor, self.acceptor, *self.cds, *self.intron]:
            if len(arr) != n:
                raise ValueError("every channel must have length n")

    @classmethod
    def zeros(cls, n):
        return cls(n)

    def favour(self, cds_ranges, intron_ranges, weight=10.0):
        """Test helper: reward CDS bases in `cds_ranges` and intronic bases in
        `intron_ranges` (half-open oriented intervals) in every phase channel;
        leave everything else at zero."""
        for a, b in cds_ranges:
            for t in range(a, b):
                for p in range(3):
                    self.cds[p][t] = weight
        for a, b in intron_ranges:
            for t in range(a, b):
                for p in range(3):
                    self.intron[p][t] = weight
        # anchor the junctions: the only rewarded donor is at each intron's
        # first base, the only rewarded acceptor at the base after its last
        if intron_ranges:
            self.donor = [-weight] * self.n
            self.acceptor = [-weight] * self.n
            for a, b in intron_ranges:
                self.donor[a] = weight
                if b < self.n:
                    self.acceptor[b] = weight
        return self
