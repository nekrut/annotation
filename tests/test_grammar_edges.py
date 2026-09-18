"""Regressions for the edge, mask and input-validation findings on PR #31:
masked intron emissions poisoning the delayed rolling sums (engels-0037),
empty input and zero-length residual introns (stalin-0044), and NaN / +inf
score channels."""
import itertools
import unittest
from math import inf, nan, log, isinf, isnan

from model.grammar import (ReferenceDecoder, DelayedEntryDecoder, DurationMixture,
                           EdgePrior, Scores, TABLES)

NEG = -inf
DECODERS = (ReferenceDecoder, DelayedEntryDecoder)


def pinned(x, start, donor, acceptor):
    """Scores that forbid U everywhere and pin one start, one donor and one
    acceptor, so the only legal chain is the one those junctions describe."""
    sc = Scores.zeros(len(x))
    for arr in (sc.u, sc.start, sc.donor, sc.acceptor):
        arr[:] = [NEG] * len(x)
    sc.start[start] = sc.donor[donor] = sc.acceptor[acceptor] = 0.0
    return sc


class MaskedIntronEmissions(unittest.TestCase):
    def test_mask_before_intron_recovers(self):
        # engels-0037: a -inf intron emission before a later legal intron
        for m in (1, 2, 20):
            x = "ATG" + "C" * m + "AAATAA"
            sc = pinned(x, 0, 3, 3 + m)
            sc.intron[0][0] = NEG                     # outside the real intron [3, 3+m)
            for cls in DECODERS:
                dec = cls(TABLES[1], DurationMixture(m=m))
                best, chains = dec.viterbi(x, sc)
                z = dec.partition(x, sc)
                self.assertAlmostEqual(z, log(0.5), places=12, msg=(m, cls))
                self.assertAlmostEqual(best, log(0.5), places=12, msg=(m, cls))
                self.assertEqual(len(chains), 1)
                c = chains[0]
                self.assertFalse(c.partial_5 or c.partial_3)
                self.assertEqual([(s.start, s.end) for s in c.introns()], [(3, 3 + m)])

    def test_mask_sweep_matches_closed_form(self):
        # engels-0037's 72 pinned-chain configurations: phase p, m, intron
        # length m or m+1, mask none/before/inside/after
        pi, q = (0.2, 0.8), (0.2, 0.6)
        dur_by_m = {m: DurationMixture(m=m, pi=(pi,) * 3, q=(q,) * 3) for m in (1, 2, 20)}
        for p, m, extra, mask in itertools.product((0, 1, 2), (1, 2, 20), (0, 1), ("none", "before", "inside", "after")):
            length = m + extra
            cut = 3 + p
            left, right = "ATGAAATAA"[:cut], "ATGAAATAA"[cut:]
            x = left + "C" * length + right
            sc = pinned(x, 0, cut, cut + length)
            for t in range(cut, cut + length):
                for ph in range(3):
                    sc.cds[ph][t] = NEG
            at = {"none": None, "before": 0, "inside": cut, "after": cut + length}[mask]
            if at is not None:
                sc.intron[p][at] = NEG
            comps = [log(pi[r]) + log(1 - q[r]) + (length - m) * log(q[r]) for r in range(2)]
            for cls in DECODERS:
                dec = cls(TABLES[1], dur_by_m[m])
                z = dec.partition(x, sc)
                best, chains = dec.viterbi(x, sc)
                msg = (p, m, length, mask, cls.__name__)
                if mask == "inside":
                    self.assertTrue(isinf(z) and z < 0, msg)
                    self.assertTrue(isinf(best) and best < 0, msg)
                    self.assertEqual(chains, [], msg)
                    continue
                self.assertAlmostEqual(z, log(sum(map(lambda v: 2.718281828459045 ** v, comps))), places=10, msg=msg)
                self.assertAlmostEqual(best, max(comps), places=12, msg=msg)
                self.assertEqual(len(chains), 1, msg)
                c = chains[0]
                self.assertEqual(c.spliced(x), "ATGAAATAA", msg)
                self.assertEqual([s.gff_phase for s in c.cds()], [0, (3 - p) % 3], msg)
                self.assertEqual([(s.start, s.end, s.component) for s in c.introns()],
                                 [(cut, cut + length, comps.index(max(comps)))], msg)

    def test_masks_leaving_window_before_a_later_donor_random(self):
        # random lattices with sparse -inf intron emissions: the delayed and
        # reference decoders must agree (finite or both -inf, never NaN)
        import random
        rng = random.Random(20260915)
        agreements = 0
        for trial in range(60):
            m, R = rng.choice([1, 2, 3, 4]), rng.choice([1, 2])
            pi = tuple(tuple(v / sum(row) for v in row) for row in [[rng.uniform(0.1, 1.0) for _ in range(R)] for _ in range(3)])
            qq = tuple(tuple(rng.uniform(0.1, 0.9) for _ in range(R)) for _ in range(3))
            dur = DurationMixture(m=m, pi=pi, q=qq)
            edges = rng.choice([None, EdgePrior(log(0.3), log(0.4))])
            n = rng.randint(2, 12)
            x = "".join(rng.choice("ACGT") for _ in range(n))
            sc = Scores.zeros(n)
            for name, arr in sc.channels():
                for t in range(n):
                    arr[t] = rng.uniform(-1.0, 1.0) + (1.0 if name.startswith("cds") else 0.0)
            for p in range(3):
                for t in range(n):
                    if rng.random() < 0.15:
                        sc.intron[p][t] = NEG
            ref, dl = ReferenceDecoder(TABLES[1], dur, edges), DelayedEntryDecoder(TABLES[1], dur, edges)
            zr, zd = ref.partition(x, sc), dl.partition(x, sc)
            (vr, cr), (vd, cd) = ref.viterbi(x, sc), dl.viterbi(x, sc)
            for a, b in ((zr, zd), (vr, vd)):
                self.assertFalse(isnan(a) or isnan(b))
                if isinf(a) or isinf(b):
                    self.assertEqual(a, b)
                else:
                    self.assertAlmostEqual(a, b, places=10)
            self.assertEqual(cr, cd)
            agreements += bool(cr)
        self.assertGreater(agreements, 10)


class EdgeEntryPolicy(unittest.TestCase):
    def test_empty_input_has_no_chain_and_unit_partition(self):
        # stalin-0044: E0 must not be terminal, so x = "" admits only the empty U path
        for cls in DECODERS:
            for edges in (EdgePrior(), EdgePrior(2.0, 2.0)):
                dec = cls(edges=edges)
                self.assertEqual(dec.partition("", Scores.zeros(0)), 0.0)
                self.assertEqual(dec.viterbi("", Scores.zeros(0)), (0.0, []))

    def test_residual_intron_needs_one_observed_base(self):
        # stalin-0044: J may not close at t = 0, so no zero-length intron is
        # ever emitted and the CDS-from-edge hypothesis is priced once (E0)
        x = "AAATAA"
        sc = Scores.zeros(len(x)).favour([(0, len(x))], [])
        sc.acceptor[0] = 10.0
        edges = EdgePrior()
        for cls in DECODERS:
            best, chains = cls(edges=edges).viterbi(x, sc)
            self.assertEqual(len(chains), 1)
            c = chains[0]
            self.assertTrue(c.partial_5 and not c.partial_3)
            self.assertEqual(c.introns(), [])
            self.assertEqual([(s.start, s.end, s.gff_phase) for s in c.cds()], [(0, 6, 0)])
            # E0("") with the phase/prefix prior, 6 favoured bases, stop at the end
            self.assertAlmostEqual(best, edges.entry - log(3.0) + 60.0, places=12)
        # with one observed intronic base the residual intron is available again
        y = "C" + x
        sc = Scores.zeros(len(y)).favour([(1, len(y))], [])
        sc.acceptor[1] = 10.0
        for cls in DECODERS:
            _, chains = cls(edges=edges).viterbi(y, sc)
            self.assertEqual([(s.start, s.end) for s in chains[0].introns()], [(0, 1)])

    def test_no_zero_length_intron_over_random_lattices(self):
        import random
        rng = random.Random(44)
        for trial in range(80):
            n = rng.randint(1, 10)
            x = "".join(rng.choice("ACGTN") for _ in range(n))
            sc = Scores.zeros(n)
            for _, arr in sc.channels():
                for t in range(n):
                    arr[t] = rng.uniform(-2.0, 2.0)
            for cls in DECODERS:
                dec = cls(TABLES[rng.choice([1, 6])], DurationMixture(m=rng.choice([1, 2, 3])), EdgePrior(1.0, 1.0))
                _, chains = dec.viterbi(x, sc)
                for c in chains:
                    for s in c.segments:
                        self.assertLess(s.start, s.end, (x, s))


class InputValidation(unittest.TestCase):
    def test_nan_and_positive_infinity_are_rejected(self):
        for cls in DECODERS:
            for bad in (nan, inf):
                for name in ("u", "start", "donor"):
                    sc = Scores.zeros(4)
                    getattr(sc, name)[2] = bad
                    with self.assertRaises(ValueError):
                        cls().partition("ATGA", sc)
                sc = Scores.zeros(4)
                sc.intron[1][0] = bad
                with self.assertRaises(ValueError):
                    cls().viterbi("ATGA", sc)

    def test_length_mismatch_is_rejected(self):
        for cls in DECODERS:
            with self.assertRaises(ValueError):
                cls().partition("ATG", Scores.zeros(4))


if __name__ == "__main__":
    unittest.main()
