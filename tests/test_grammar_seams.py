"""Proposal 3.3 / 3.4: an intron spanning a checkpoint or chunk seam keeps
its prefix, duration component and score, and a seam offers no new
start/end permission. The delayed decoder is checkpointed at every seam
and resumed from the checkpoint alone; results must equal the unchunked
run and the expanded reference decoder."""
import itertools
import random
import unittest
from math import log, inf, isinf

from model.grammar import (ReferenceDecoder, DelayedEntryDecoder, DurationMixture,
                           EdgePrior, Scores, TABLES)

M = 20
INTRON = "GT" + "A" * (M - 4) + "AG"


def same(test, a, b, msg=None):
    if isinf(a) or isinf(b):
        test.assertEqual(a, b, msg)
    else:
        test.assertAlmostEqual(a, b, places=10, msg=msg)


class Seams(unittest.TestCase):
    def test_intron_spanning_every_seam_position(self):
        # the phase-1 acceptance case; seams at every interior boundary,
        # including inside the mandatory interval and at the donor/acceptor
        x = "ATGAAAT" + INTRON + "AA"
        n = len(x)
        sc = Scores.zeros(n).favour([(0, 7), (7 + M, n)], [(7, 7 + M)])
        for edges in (None, EdgePrior(log(0.3), log(0.4))):
            dec = DelayedEntryDecoder(TABLES[1], DurationMixture(m=M, pi=((0.7, 0.3),) * 3, q=((0.9, 0.5),) * 3), edges)
            z0 = dec.partition(x, sc)
            v0, c0 = dec.viterbi(x, sc)
            self.assertEqual(c0[0].spliced(x), "ATGAAATAA")
            for t in range(1, n):
                same(self, dec.partition(x, sc, seams=[t]), z0, t)
                v, c = dec.viterbi(x, sc, seams=[t])
                same(self, v, v0, t)
                self.assertEqual(c, c0, t)
            # many seams at once
            same(self, dec.partition(x, sc, seams=range(1, n)), z0)
            v, c = dec.viterbi(x, sc, seams=range(1, n))
            same(self, v, v0)
            self.assertEqual(c, c0)

    def test_seam_grants_no_partial_entry_or_exit(self):
        # with edges on, a CDS that is only favoured after the seam must still
        # be entered through U/S (or E0 at boundary 0), never at the seam
        x = "CCCCC" + "ATGAAATAA"
        n = len(x)
        sc = Scores.zeros(n).favour([(5, n)], [])
        for t in range(5):
            sc.u[t] = -inf                  # U forbidden before the CDS: only E0/J entry at 0 remains
            for p in range(3):              # break ties among the equal-score edge hypotheses
                sc.cds[p][t] = -0.1 * (t + 1) * (p + 1)
        dec = DelayedEntryDecoder(TABLES[1], DurationMixture(m=2), EdgePrior(log(0.3), log(0.4)))
        ref = ReferenceDecoder(TABLES[1], DurationMixture(m=2), EdgePrior(log(0.3), log(0.4)))
        vr, cr = ref.viterbi(x, sc)
        for seams in ([5], [4, 5, 6], list(range(1, n))):
            v, c = dec.viterbi(x, sc, seams=seams)
            same(self, v, vr, seams)
            self.assertEqual(c, cr, seams)
            self.assertTrue(c[0].partial_5)  # entered at boundary 0, not at the seam
            same(self, dec.partition(x, sc, seams=seams), ref.partition(x, sc), seams)

    def test_random_lattices_with_random_seams_match_reference(self):
        rng = random.Random(2026091518)
        chains = 0
        for trial in range(150):
            m, R = rng.choice([1, 2, 3, 4]), rng.choice([1, 2, 3])
            pi = tuple(tuple(v / sum(row) for v in row) for row in [[rng.uniform(0.1, 1.0) for _ in range(R)] for _ in range(3)])
            q = tuple(tuple(rng.uniform(0.1, 0.9) for _ in range(R)) for _ in range(3))
            dur = DurationMixture(m=m, pi=pi, q=q)
            edges = rng.choice([None, EdgePrior(log(0.3), log(0.4))])
            table = rng.choice([1, 6])
            n = rng.randint(2, 12)
            x = "".join(rng.choice("ACGTACGTNRY") for _ in range(n))
            sc = Scores.zeros(n)
            for name, arr in sc.channels():
                for t in range(n):
                    arr[t] = rng.uniform(-1.0, 1.0) + (1.0 if name.startswith("cds") else 0.0)
                    if name.startswith("intron") and rng.random() < 0.1:
                        arr[t] = -inf
            seams = sorted(t for t in range(1, n) if rng.random() < 0.4)
            ref = ReferenceDecoder(TABLES[table], dur, edges)
            dl = DelayedEntryDecoder(TABLES[table], dur, edges)
            zr, vr_c = ref.partition(x, sc), ref.viterbi(x, sc)
            same(self, dl.partition(x, sc, seams=seams), zr, (x, seams))
            vd, cd = dl.viterbi(x, sc, seams=seams)
            same(self, vd, vr_c[0], (x, seams))
            self.assertEqual(cd, vr_c[1], (x, seams))
            chains += bool(cd)
        self.assertGreater(chains, 30)

    def test_seams_must_be_interior(self):
        dec = DelayedEntryDecoder()
        for bad in ([0], [3], [-1]):
            with self.assertRaises(ValueError):
                dec.partition("ATG", Scores.zeros(3), seams=bad)


if __name__ == "__main__":
    unittest.main()
