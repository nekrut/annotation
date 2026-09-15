"""Reverse-complement coordinate and phase agreement (proposal 3.1: GFF3
phase derived along the traced chain, with strand-aware coordinates)."""
import unittest

from model.grammar import (ReferenceDecoder, DelayedEntryDecoder, DurationMixture, Scores,
                           TABLES, reverse_complement, genomic_features, gff3_rows)

M = 20
INTRON = "GT" + "A" * (M - 4) + "AG"


class ReverseComplement(unittest.TestCase):
    def test_iupac_round_trip(self):
        x = "ACGTRYSWKMBDHVNX"
        self.assertEqual(reverse_complement(reverse_complement(x)), x)
        self.assertEqual(reverse_complement("ATGAAATAA"), "TTATTTCAT")
        self.assertEqual(reverse_complement("RYKMBDHV"), "BDHVKMRY")

    def decode_both(self, oriented, sc):
        out = []
        for cls in (ReferenceDecoder, DelayedEntryDecoder):
            score, chains = cls(TABLES[1], DurationMixture(m=M)).viterbi(oriented, sc)
            self.assertEqual(len(chains), 1)
            out.append((score, chains[0]))
        self.assertAlmostEqual(out[0][0], out[1][0], places=12)
        self.assertEqual(out[0][1], out[1][1])
        return out[0]

    def test_minus_strand_gene_maps_back_with_phase(self):
        # the phase-1 acceptance case planted on the minus strand of a window
        oriented = "ATGAAAT" + INTRON + "AA"
        n = len(oriented)
        genome = reverse_complement(oriented)
        self.assertEqual(reverse_complement(genome), oriented)
        cds_o, intron_o = [(0, 7), (7 + M, n)], [(7, 7 + M)]
        sc = Scores.zeros(n).favour(cds_o, intron_o)
        score, chain = self.decode_both(reverse_complement(genome), sc)
        self.assertEqual(chain.spliced(oriented), "ATGAAATAA")
        feats = genomic_features(chain, n, "-")
        self.assertEqual([(f.kind, f.start, f.end, f.phase) for f in feats],
                         [("cds", 0, 2, 2), ("intron", 2, 2 + M, None), ("cds", 2 + M, n, 0)])
        # the split stop's AA sits at the genomic low end as TT (its T is
        # across the intron); the initiator ATG at the high end as CAT
        self.assertEqual(genome[0:2], "TT")
        self.assertEqual(genome[n - 3:n], "CAT")
        # the same chain on the plus strand is the identity map
        plus = genomic_features(chain, n, "+")
        self.assertEqual([(f.start, f.end, f.phase) for f in plus if f.kind == "cds"], [(0, 7, 0), (7 + M, n, 2)])
        rows = gff3_rows(chain, "chrT", n, "-", "g1")
        self.assertEqual(rows, ["chrT\tgrammar\tCDS\t1\t2\t.\t-\t2\tParent=g1",
                                f"chrT\tgrammar\tCDS\t{3 + M}\t{n}\t.\t-\t0\tParent=g1"])
        # spliced genomic minus-strand CDS re-reads as the oriented CDS
        spliced = "".join(genome[f.start:f.end] for f in feats if f.kind == "cds")
        self.assertEqual(reverse_complement(spliced), "ATGAAATAA")

    def test_gff_phase_consistent_with_spliced_length_on_both_strands(self):
        # phase-2 case: ATGAAATA |intron| A ; GFF phase of the second CDS is 1
        oriented = "ATGAAATA" + INTRON + "A"
        n = len(oriented)
        sc = Scores.zeros(n).favour([(0, 8), (8 + M, n)], [(8, 8 + M)])
        _, chain = self.decode_both(oriented, sc)
        for strand in "+-":
            cds = [f for f in genomic_features(chain, n, strand) if f.kind == "cds"]
            # walking the chain 5' to 3' on either strand, phase_k = (3 - emitted) % 3
            ordered = cds if strand == "+" else list(reversed(cds))
            emitted = 0
            for f in ordered:
                self.assertEqual(f.phase, (3 - emitted % 3) % 3)
                emitted += f.end - f.start
            self.assertEqual(emitted % 3, 0)


if __name__ == "__main__":
    unittest.main()
