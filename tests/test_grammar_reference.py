"""Proposal section 3.4 acceptance cases against the expanded reference decoder."""
import unittest
from math import log, isinf

from model.grammar import ReferenceDecoder, DurationMixture, Scores, TABLES

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
    """Every legal path from U at 0 to U at n with its score (brute force)."""
    from model.grammar.reference import U
    out = []

    def dfs(state, t, score, path):
        if t == len(x):
            if state == U:
                out.append((score, tuple(path)))
            return
        for nxt, s, _ in dec.transitions(state, x, t, sc):
            if isinf(s):
                continue
            dfs(nxt, t + 1, score + s, path + [nxt])
    dfs(U, 0, 0.0, [U])
    return out


class BruteForceAgreement(unittest.TestCase):
    """Section 3.4: brute-force agreement of the CRF partition and of Viterbi."""

    def test_partition_and_viterbi_match_path_enumeration(self):
        import random
        from math import exp
        rng = random.Random(20260915)
        m = 3
        dur = DurationMixture(m=m, pi=((0.6, 0.4),) * 3, q=((0.3, 0.8),) * 3)
        dec = ReferenceDecoder(TABLES[1], dur)
        for trial in range(6):
            n = rng.randint(8, 13)
            x = "".join(rng.choice("ACGTACGTN") for _ in range(n))
            sc = Scores.zeros(n)
            for arr in [sc.u, sc.start, sc.stop, sc.donor, sc.acceptor, *sc.cds, *sc.intron]:
                for t in range(n):
                    arr[t] = rng.uniform(-2.0, 2.0)
            paths = enumerate_paths(dec, x, sc)
            self.assertTrue(paths)
            total = max(s for s, _ in paths)
            total = total + log(sum(exp(s - total) for s, _ in paths))
            self.assertAlmostEqual(dec.partition(x, sc), total, places=9, msg=x)
            best_score, _ = dec.viterbi(x, sc)
            self.assertAlmostEqual(best_score, max(s for s, _ in paths), places=9, msg=x)


if __name__ == "__main__":
    unittest.main()
