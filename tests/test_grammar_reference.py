"""Proposal section 3.4 acceptance cases against the expanded reference decoder."""
import unittest
from math import log, exp, isinf

from model.grammar import (ReferenceDecoder, DelayedEntryDecoder, DurationMixture,
                           EdgePrior, Scores, TABLES)

M = 20
INTRON = "GT" + "A" * (M - 4) + "AG"          # exactly M bases
assert len(INTRON) == M


def build(*parts):
    """Concatenate ('cds', s) / ('intron', s) parts; return sequence and
    favouring Scores so the intended path wins over all-U."""
    x, cds, introns, t = "", [], [], 0
    for kind, s in parts:
        (cds if kind == "cds" else introns).append((t, t + len(s)))
        x += s
        t += len(s)
    sc = Scores.zeros(len(x)).favour(cds, introns)
    return x, sc, cds, introns


class AcceptanceCases(unittest.TestCase):
    def decode(self, x, sc, table=1):
        dec = ReferenceDecoder(TABLES[table], DurationMixture(m=M))
        score, chains = dec.viterbi(x, sc)
        return score, chains

    def test_phase1_intron_keeps_prefix_T(self):
        x, sc, cds, introns = build(("cds", "ATGAAAT"), ("intron", INTRON), ("cds", "AA"))
        _, chains = self.decode(x, sc)
        self.assertEqual(len(chains), 1)
        c = chains[0]
        self.assertEqual(c.spliced(x), "ATGAAATAA")
        self.assertEqual([(s.start, s.end) for s in c.cds()], cds)
        self.assertEqual([(s.start, s.end) for s in c.introns()], introns)
        self.assertEqual([s.gff_phase for s in c.cds()], [0, 2])   # second exon starts after p=1

    def test_phase2_intron_keeps_prefix_TA(self):
        x, sc, cds, _ = build(("cds", "ATGAAATA"), ("intron", INTRON), ("cds", "A"))
        _, chains = self.decode(x, sc)
        self.assertEqual(chains[0].spliced(x), "ATGAAATAA")
        self.assertEqual([s.gff_phase for s in chains[0].cds()], [0, 1])

    def test_split_initiator_through_S_states(self):
        x, sc, cds, _ = build(("cds", "A"), ("intron", INTRON), ("cds", "TGAAATAA"))
        _, chains = self.decode(x, sc)
        self.assertEqual(len(chains), 1)
        self.assertEqual(chains[0].spliced(x), "ATGAAATAA")
        self.assertEqual([(s.start, s.end) for s in chains[0].cds()], cds)

    def test_table6_reads_through_TAA_table1_stops(self):
        x, sc, _, _ = build(("cds", "ATGAAAT"), ("intron", INTRON), ("cds", "AATGA"))
        _, c6 = self.decode(x, sc, table=6)
        self.assertEqual(c6[0].spliced(x), "ATGAAATAATGA")
        _, c1 = self.decode(x, sc, table=1)
        self.assertEqual(c1[0].spliced(x), "ATGAAATAA")   # terminates at the earlier TAA

    def test_intron_lengths_m_minus_1_m_m_plus_1(self):
        dur = DurationMixture(m=M, pi=((0.7, 0.3),) * 3, q=((0.9, 0.5),) * 3)
        for length in (M - 1, M, M + 1):
            intron = "G" * length
            x = "ATGAAAT" + intron + "AA"
            sc = Scores.zeros(len(x))
            # forbid U everywhere and pin the junctions, so the only legal
            # chain is ATGAAAT |intron| AA (an unpinned intron could shift one
            # base right and end on TGA instead)
            for t in range(len(x)):
                sc.u[t] = float("-inf")
                sc.donor[t] = float("-inf")
                sc.acceptor[t] = float("-inf")
            sc.donor[7] = 0.0
            sc.acceptor[7 + length] = 0.0
            dec = ReferenceDecoder(TABLES[1], dur)
            score, chains = dec.viterbi(x, sc)
            if length < M:
                self.assertTrue(isinf(score))
                continue
            p = 1   # prefix T is pending across the intron
            per_r = [dur.log_pi(p, r) + (length - M) * dur.log_q(p, r) + dur.log_1mq(p, r) for r in range(2)]
            self.assertAlmostEqual(score, max(per_r), places=12)
            self.assertEqual(chains[0].introns()[0].component, per_r.index(max(per_r)))
            self.assertAlmostEqual(dec.partition(x, sc), dur.log_prob(p, length), places=12)

    def test_no_gene_when_nothing_favoured(self):
        x = "ATGAAATAA"
        dec = ReferenceDecoder(TABLES[1], DurationMixture(m=M))
        score, chains = dec.viterbi(x, Scores.zeros(len(x)))
        self.assertEqual(score, 0.0)          # all-U path scores zero and ties are broken by first arrival
        self.assertEqual(chains, [])

    def test_ambiguous_base_marks_chain_uncertain_and_uses_prior(self):
        x, sc, _, _ = build(("cds", "ATGANATAA"),)
        dec = ReferenceDecoder(TABLES[1], DurationMixture(m=M))
        score, chains = dec.viterbi(x, sc)
        self.assertTrue(chains[0].uncertain)
        # 9 favoured CDS bases at weight 10, one uniform prior over four bases
        self.assertAlmostEqual(score, 90.0 + log(0.25), places=12)


def enumerate_paths(dec, x, sc):
    """Every legal path from boundary 0 to boundary n with its score (brute
    force over the reference transition generator, initial and terminal
    weights included, so edge partials are enumerated too)."""
    out = []

    def dfs(state, t, score, path):
        if t == len(x):
            fin = dec.terminal(state)
            if not isinf(fin):
                out.append((score + fin, tuple(path)))
            return
        for nxt, s, _ in dec.transitions(state, x, t, sc):
            if isinf(s):
                continue
            dfs(nxt, t + 1, score + s, path + [nxt])
    for state, w in dec.initial().items():
        dfs(state, 0, w, [state])
    return out


def random_scores(rng, n, cds_bonus=0.0):
    sc = Scores.zeros(n)
    for arr in [sc.u, sc.start, sc.stop, sc.donor, sc.acceptor, *sc.cds, *sc.intron]:
        for t in range(n):
            arr[t] = rng.uniform(-2.0, 2.0)
    for arr in sc.cds:
        for t in range(n):
            arr[t] += cds_bonus
    return sc


class BruteForceAgreement(unittest.TestCase):
    """Section 3.4: brute-force agreement of the CRF partition and of Viterbi."""

    def test_partition_and_viterbi_match_path_enumeration(self):
        import random
        from math import exp
        rng = random.Random(20260915)
        m = 3
        dur = DurationMixture(m=m, pi=((0.6, 0.4),) * 3, q=((0.3, 0.8),) * 3)
        for trial in range(6):
            n = rng.randint(8, 13)
            x = "".join(rng.choice("ACGTACGTN") for _ in range(n))
            sc = Scores.zeros(n)
            for arr in [sc.u, sc.start, sc.stop, sc.donor, sc.acceptor, *sc.cds, *sc.intron]:
                for t in range(n):
                    arr[t] = rng.uniform(-2.0, 2.0)
            for edges in (None, EdgePrior(log(0.3), log(0.4))):
                dec = ReferenceDecoder(TABLES[1], dur, edges)
                paths = enumerate_paths(dec, x, sc)
                self.assertTrue(paths)
                total = max(s for s, _ in paths)
                total = total + log(sum(exp(s - total) for s, _ in paths))
                self.assertAlmostEqual(dec.partition(x, sc), total, places=9, msg=x)
                best_score, _ = dec.viterbi(x, sc)
                self.assertAlmostEqual(best_score, max(s for s, _ in paths), places=9, msg=x)


class ReviewerFindings(unittest.TestCase):
    """Regressions for the interim findings on PR #31 (engels-0036, stalin-0043)."""

    def test_one_base_intron_at_m1_keeps_its_component(self):
        # engels-0036: with m = 1 the donor enters a tail directly and the
        # traceback dropped the component; both lengths must report r = 1
        dur = DurationMixture(m=1, pi=((0.2, 0.8),) * 3, q=((0.2, 0.6),) * 3)
        for length in (1, 2):
            x = "ATG" + "C" * length + "AAATAA"
            sc = Scores.zeros(len(x))
            for arr in (sc.u, sc.start, sc.donor, sc.acceptor):
                arr[:] = [float("-inf")] * len(x)
            sc.start[0] = sc.donor[3] = sc.acceptor[3 + length] = 0.0
            per_r = [dur.log_pi(0, r) + dur.log_1mq(0, r) + (length - 1) * dur.log_q(0, r) for r in range(2)]
            for cls in (ReferenceDecoder, DelayedEntryDecoder):
                dec = cls(TABLES[1], dur)
                best, chains = dec.viterbi(x, sc)
                self.assertAlmostEqual(best, max(per_r), places=12)
                self.assertAlmostEqual(dec.partition(x, sc), dur.log_prob(0, length), places=12)
                intron = chains[0].introns()[0]
                self.assertEqual((intron.start, intron.end), (3, 3 + length))
                self.assertEqual(intron.component, per_r.index(max(per_r)), (cls.__name__, length))

    def test_duration_mixture_rejects_bad_shapes_and_weights(self):
        # stalin-0043: unequal component counts silently dropped mass; zero or
        # negative weights constructed and then failed in log_pi
        with self.assertRaises(ValueError):
            DurationMixture(m=2, pi=((1.0,), (0.5, 0.5), (1.0,)), q=((0.5,), (0.5, 0.5), (0.5,)))
        with self.assertRaises(ValueError):
            DurationMixture(m=2, pi=((1.0, 0.0),) * 3, q=((0.5, 0.5),) * 3)
        with self.assertRaises(ValueError):
            DurationMixture(m=2, pi=((1.2, -0.2),) * 3, q=((0.5, 0.5),) * 3)
        with self.assertRaises(ValueError):
            DurationMixture(m=2, pi=((0.5, 0.5),) * 2, q=((0.5, 0.5),) * 2)
        with self.assertRaises(ValueError):
            DurationMixture(m=2, pi=((0.5, 0.5),) * 3, q=((0.5, 1.0),) * 3)
        d = DurationMixture(m=2, pi=((0.5, 0.5),) * 3, q=((0.5, 0.5),) * 3)
        self.assertAlmostEqual(sum(exp(d.log_prob(1, l)) for l in range(2, 200)), 1.0, places=12)


class EdgePartials(unittest.TestCase):
    """Proposal 3.1: E/I entry and E/S/I exit at an actual sequence edge with
    an explicit partial-end prior; paths with no observed CDS are discarded."""
    EDGES = EdgePrior(log(0.5), log(0.5))

    def decoders(self, table=1, m=M):
        dur = DurationMixture(m=m)
        return [cls(TABLES[table], dur, self.EDGES) for cls in (ReferenceDecoder, DelayedEntryDecoder)]

    def test_partial_5_cds_from_the_edge_with_prefix_prior(self):
        # sequence starts inside a codon: (AA)A AAA TAA, two bases pending
        # upstream; the stop is rewarded so that completing the chain beats
        # the cheaper censored exit E0('') AAA ATA A|
        x = "AAAATAA"
        sc = Scores.zeros(len(x)).favour([(0, len(x))], [])
        sc.stop[6] = 10.0
        for dec in self.decoders():
            score, chains = dec.viterbi(x, sc)
            c = chains[0]
            self.assertTrue(c.partial_5)
            self.assertFalse(c.partial_3)
            self.assertEqual([(s.start, s.end, s.gff_phase) for s in c.cds()], [(0, 7, 1)])
            # 7 favoured bases, the stop, entry prior, and the p = 2 prefix prior 1/3 * 1/16
            self.assertAlmostEqual(score, 80.0 + log(0.5) - log(3.0) - 2 * log(4.0), places=12)
        # without the stop reward the censored exit wins, with the p = 0 prior
        for dec in self.decoders():
            score, chains = dec.viterbi(x, Scores.zeros(len(x)).favour([(0, len(x))], []))
            self.assertTrue(chains[0].partial_3)
            self.assertAlmostEqual(score, 70.0 + 2 * log(0.5) - log(3.0), places=12)

    def test_partial_3_exit_inside_an_intron_and_inside_a_codon(self):
        x, sc, cds, introns = build(("cds", "ATGAAAT"), ("intron", INTRON))
        for dec in self.decoders():
            _, chains = dec.viterbi(x, sc)
            c = chains[0]
            self.assertTrue(c.partial_3)
            self.assertEqual([(s.start, s.end) for s in c.introns()], introns)
            self.assertEqual(c.introns()[0].prefix_len, 1)        # prefix T is pending
        # censored intron shorter than m: the delayed decoder must offer the
        # pending-donor exit that the reference reaches through I(c, k)
        x, sc, cds, introns = build(("cds", "ATGAAAT"), ("intron", INTRON[: M - 5]))
        for dec in self.decoders():
            _, chains = dec.viterbi(x, sc)
            self.assertTrue(chains[0].partial_3)
            self.assertEqual([(s.start, s.end) for s in chains[0].introns()], introns)
            self.assertIsNone(chains[0].introns()[0].component)
        x, sc, _, _ = build(("cds", "ATGAAAT"),)
        for dec in self.decoders():
            _, chains = dec.viterbi(x, sc)
            self.assertTrue(chains[0].partial_3)
            self.assertEqual(chains[0].spliced(x), "ATGAAAT")

    def test_residual_intron_entry_then_exit_into_cds(self):
        # sequence starts inside an intron; first observed CDS base closes it
        # and (A)AA TAA completes with the rewarded stop
        x, sc, cds, introns = build(("intron", INTRON), ("cds", "AATAA"))
        sc.stop[len(x) - 1] = 10.0
        for dec in self.decoders():
            _, chains = dec.viterbi(x, sc)
            c = chains[0]
            self.assertTrue(c.partial_5)
            self.assertEqual([(s.start, s.end) for s in c.introns()], introns)
            self.assertEqual(c.spliced(x), "AATAA")
            self.assertEqual(c.introns()[0].prefix_len, 1)      # E('A') held across the residual intron

    def test_no_observed_cds_is_discarded(self):
        # CDS forbidden everywhere: the only chains would be residual intron
        # in, residual intron out, which is not a chain
        x = INTRON
        sc = Scores.zeros(len(x)).favour([], [(0, len(x))])
        for p in range(3):
            sc.cds[p] = [float("-inf")] * len(x)
        for dec in self.decoders():
            score, chains = dec.viterbi(x, sc)
            self.assertEqual(chains, [])
            self.assertEqual(score, 0.0)      # the all-U path
        # no edges: the same favouring cannot produce a partial either
        x = "AAAATAA"
        sc = Scores.zeros(len(x)).favour([(0, len(x))], [])
        _, chains = ReferenceDecoder(TABLES[1], DurationMixture(m=M)).viterbi(x, sc)
        self.assertEqual(chains, [])


class DelayedEntryAgreement(unittest.TestCase):
    """Proposal 3.2/3.4: forward and Viterbi agreement between the expanded
    reference and the delayed-entry recurrence, including the traced chains."""

    def check(self, table, dur, edges, x, sc):
        ref = ReferenceDecoder(TABLES[table], dur, edges)
        dl = DelayedEntryDecoder(TABLES[table], dur, edges)
        pr, pd = ref.partition(x, sc), dl.partition(x, sc)
        vr, cr = ref.viterbi(x, sc)
        vd, cd = dl.viterbi(x, sc)
        msg = (table, dur.m, dur.R, edges is not None, x)
        if isinf(pr) or isinf(pd):
            self.assertEqual(pr, pd, msg)
        else:
            self.assertAlmostEqual(pr, pd, places=10, msg=msg)
        if isinf(vr) or isinf(vd):
            self.assertEqual(vr, vd, msg)
        else:
            self.assertAlmostEqual(vr, vd, places=10, msg=msg)
        self.assertEqual(cr, cd, msg)
        return bool(cr)

    def test_exhaustive_tiny_lattices(self):
        # every sequence over ACGT of length 1..5, m = 2, R = 2, edges on,
        # one fixed random score set per length
        import itertools, random
        rng = random.Random(2026091513)
        dur = DurationMixture(m=2, pi=((0.6, 0.4), (0.3, 0.7), (0.5, 0.5)), q=((0.3, 0.8), (0.5, 0.6), (0.2, 0.9)))
        edges = EdgePrior(log(0.3), log(0.4))
        chains = 0
        for n in range(1, 6):
            sc = random_scores(rng, n, cds_bonus=1.0)
            for bases in itertools.product("ACGT", repeat=n):
                chains += self.check(1, dur, edges, "".join(bases), sc)
        self.assertGreater(chains, 0)

    def test_random_lattices_with_ambiguity_tables_and_minimum_lengths(self):
        import random
        rng = random.Random(20260915)
        chains = 0
        for trial in range(120):
            m, R = rng.choice([1, 2, 3, 4]), rng.choice([1, 2, 3])
            pi = tuple(tuple(v / sum(row) for v in row)
                       for row in [[rng.uniform(0.1, 1.0) for _ in range(R)] for _ in range(3)])
            q = tuple(tuple(rng.uniform(0.1, 0.9) for _ in range(R)) for _ in range(3))
            dur = DurationMixture(m=m, pi=pi, q=q)
            edges = rng.choice([None, EdgePrior(log(0.3), log(0.4))])
            n = rng.randint(1, 12)
            x = "".join(rng.choice("ACGTACGTNRY") for _ in range(n))
            chains += self.check(rng.choice([1, 6]), dur, edges, x, random_scores(rng, n, cds_bonus=1.0))
        self.assertGreater(chains, 20)


if __name__ == "__main__":
    unittest.main()
