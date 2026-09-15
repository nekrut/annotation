"""Proposal 3.6 adapter fixtures: a complete split codon, a nonzero-phase
edge partial, an interior partial, an internal range tag, a short-gap
frame adjustment, a translation exception, an N, and a declaration/
sequence disagreement, on both strands, plus the representative and
topology policies of 4.4. Each fixture is planted in a synthetic genome
and audited through the same code path as the train panel."""
import gzip
import os
import tempfile
import unittest

from model.labels.admission import (audit_species, metadata_flags, sequence_flags, revcomp,
                                    _score_module, load_gff_rows)

MIN_LEN = _score_module().MIN_SEQ_LEN


class Genome:
    """Builds a two-sequence genome (names that the scorer keeps: chrM would be dropped as an organelle) and its GFF3 from planted CDS pieces."""

    def __init__(self):
        self.seqs = {"chrP": ["C"] * max(MIN_LEN, 5000), "chrB": ["C"] * max(MIN_LEN, 5000)}
        self.lines = []
        self.n = 0

    def gene(self, seqid, strand, pos, pieces, gaps, phases=None, attrs="", tx_attrs="",
             biotype="protein_coding", gene_id=None, pseudo=False, tid=None, row_attrs=None):
        """Plant CDS `pieces` (transcriptional order) separated by `gaps`
        genomic bases from oriented coordinate `pos`; phases follow the
        spliced frame unless given. `attrs` go on the first row in
        transcriptional order, `row_attrs` (a list) on each row. Returns the
        transcript id, which `tid` may force (to plant a repeated id)."""
        self.n += 1
        tid, gid = tid or "t%d" % self.n, gene_id or "g%d" % self.n
        L = len(self.seqs[seqid])
        # oriented coordinates: 0-based along the strand
        cur = pos
        rows = []
        p = 0
        for i, piece in enumerate(pieces):
            rows.append((cur, cur + len(piece), piece, (-p) % 3 if phases is None else phases[i]))
            p = (p + len(piece)) % 3
            cur += len(piece) + (gaps[i] if i < len(gaps) else 0)
        genomic = []
        for s, e, piece, ph in rows:
            if strand == "+":
                gs, ge = s, e
                self.seqs[seqid][gs:ge] = list(piece)
            else:
                gs, ge = L - e, L - s
                self.seqs[seqid][gs:ge] = list(revcomp(piece))
            genomic.append((gs + 1, ge, ph))
        span = min(g[0] for g in genomic), max(g[1] for g in genomic)
        ftype = "pseudogene" if pseudo else "gene"
        self.lines.append("%s\tt\t%s\t%d\t%d\t.\t%s\t.\tID=%s;gene_biotype=%s" % (seqid, ftype, span[0], span[1], strand, gid, biotype))
        self.lines.append("%s\tt\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s%s" % (seqid, span[0], span[1], strand, tid, gid, tx_attrs))
        for i, (gs, ge, ph) in enumerate(genomic):
            # row attributes go on the first row in transcriptional order only
            self.lines.append("%s\tt\texon\t%d\t%d\t.\t%s\t.\tParent=%s" % (seqid, gs, ge, strand, tid))
            extra = (attrs if i == 0 else "") + (row_attrs[i] if row_attrs else "")
            self.lines.append("%s\tt\tCDS\t%d\t%d\t.\t%s\t%s\tID=cds-%s;Parent=%s%s" % (seqid, gs, ge, strand, ph, tid, tid, extra))
        return tid

    def write(self, d, order=("chrP", "chrB")):
        gff = os.path.join(d, "x.gff")
        fna = os.path.join(d, "x.fna.gz")
        with open(gff, "w") as fh:
            fh.write("##gff-version 3\n")
            for name, s in self.seqs.items():
                fh.write("##sequence-region %s 1 %d\n" % (name, len(s)))
            fh.write("\n".join(self.lines) + "\n")
        with gzip.open(fna, "wt") as fh:
            for name in order:
                fh.write(">%s\n%s\n" % (name, "".join(self.seqs[name])))
        return gff, fna


def five_partial(strand):
    return ";start_range=.,1" if strand == "+" else ";end_range=1,."


def three_partial(strand):
    return ";end_range=1,." if strand == "+" else ";start_range=.,1"


class Admission(unittest.TestCase):
    def audit(self, g, order=("chrP", "chrB")):
        with tempfile.TemporaryDirectory() as d:
            gff, fna = g.write(d, order)
            res = audit_species("fixture", gff, fna, m=20, table=1)
        return {r["transcript"]: r for r in res.manifest_rows}, res.summary

    def audit_keyed(self, g, order=("chrP", "chrB")):
        """Rows keyed on the full identity, for repeated transcript ids."""
        with tempfile.TemporaryDirectory() as d:
            gff, fna = g.write(d, order)
            res = audit_species("fixture", gff, fna, m=20, table=1)
        return {(r["seqid"], r["strand"], r["transcript"]): r for r in res.manifest_rows}, res.summary

    def test_fixtures_on_both_strands(self):
        for strand, seqid in (("+", "chrP"), ("-", "chrB")):
            g = Genome()
            pos = 100
            ids = {}
            # complete split codon across a legal intron: ATG AAA T|intron|AA
            ids["split"] = g.gene(seqid, strand, pos, ["ATGAAAT", "AA"], [25]); pos += 200
            # nonzero-phase 5' partial at the true sequence edge: two bases of a
            # codon observed, then in frame to a stop
            ids["edge"] = g.gene(seqid, strand, 0, ["AAGCCTAA"], [], phases=[2],
                                 attrs=";start_range=.,1" if strand == "+" else ";end_range=1,.")
            # the same declaration away from the edge is masked
            ids["edge_away"] = g.gene(seqid, strand, pos, ["AAGCCTAA"], [], phases=[2],
                                      attrs=";start_range=.,1" if strand == "+" else ";end_range=1,."); pos += 200
            # interior partial: partial=true without a range
            ids["interior"] = g.gene(seqid, strand, pos, ["ATGAAA", "TAA"], [30], attrs=";partial=true"); pos += 200
            # internal range tag on a middle row
            ids["internal_range"] = g.gene(seqid, strand, pos, ["ATG", "AAA", "TAA"], [30, 30]); pos += 200
            g.lines[-3] += ";start_range=.,%d" % 1
            # short-gap frame adjustment: one skipped base between rows
            ids["short_gap"] = g.gene(seqid, strand, pos, ["ATGAAA", "GCCTAA"], [1]); pos += 200
            # translation exception and a CDS exception tag
            ids["transl"] = g.gene(seqid, strand, pos, ["ATGAAATAA"], [], attrs=";transl_except=(pos:1..3,aa:Met)"); pos += 200
            ids["exc"] = g.gene(seqid, strand, pos, ["ATGAAATAA"], [], attrs=";exception=ribosomal slippage"); pos += 200
            # an N inside the CDS
            ids["n"] = g.gene(seqid, strand, pos, ["ATGANATAA"], []); pos += 200
            # declaration/sequence disagreements: no initiator, internal stop, no stop, bad length
            ids["no_init"] = g.gene(seqid, strand, pos, ["CTGAAATAA"], []); pos += 200
            ids["internal_stop"] = g.gene(seqid, strand, pos, ["ATGTAAAAATAA"], []); pos += 200
            ids["no_stop"] = g.gene(seqid, strand, pos, ["ATGAAAAAA"], []); pos += 200
            ids["frame"] = g.gene(seqid, strand, pos, ["ATGAAAATAA"], [], phases=[0]); pos += 200
            # inconsistent phase on the second row
            ids["phase"] = g.gene(seqid, strand, pos, ["ATGAAAT", "AA"], [25], phases=[0, 0]); pos += 200
            # overlapping rows
            ids["overlap"] = g.gene(seqid, strand, pos, ["ATGAAA", "AAATAA"], [-2]); pos += 200
            # a pseudogene and a non-coding biotype are filtered before selection
            ids["pseudo"] = g.gene(seqid, strand, pos, ["ATGAAATAA"], [], pseudo=True, biotype="pseudogene"); pos += 200
            ids["nc"] = g.gene(seqid, strand, pos, ["ATGAAATAA"], [], biotype="lncRNA"); pos += 200
            # two isoforms of one gene: the longer CDS represents
            ids["iso_short"] = g.gene(seqid, strand, pos, ["ATGAAATAA"], [], gene_id="shared"); 
            ids["iso_long"] = g.gene(seqid, strand, pos + 300, ["ATGAAAAAATAA"], [], gene_id="shared"); pos += 700
            # two genes whose spans overlap: both masked by topology
            ids["ov1"] = g.gene(seqid, strand, pos, ["ATGAAA", "TAA"], [60])
            ids["ov2"] = g.gene(seqid, strand, pos + 20, ["ATGCCCTAA"], []); pos += 300
            rows, summary = self.audit(g)
            r = lambda k: rows[ids[k]]
            self.assertEqual(r("split")["status"], "admitted", strand)
            self.assertEqual(r("split")["reasons"], "-")
            self.assertEqual(r("edge")["status"], "admitted", strand)
            self.assertEqual(r("edge")["missing_prefix"], 1)
            self.assertIn("partial_5", r("edge")["reasons"].split(","))
            self.assertEqual(r("edge_away")["status"], "masked")
            self.assertIn("partial_away_from_edge", r("edge_away")["reasons"])
            self.assertIn("partial_unlocated", r("interior")["reasons"])
            self.assertIn("partial_unlocated", r("internal_range")["reasons"])
            self.assertIn("short_gap", r("short_gap")["reasons"])
            self.assertNotIn("phase_inconsistent", r("short_gap")["reasons"])
            self.assertIn("transl_except", r("transl")["reasons"])
            self.assertIn("exception", r("exc")["reasons"])
            self.assertEqual(r("n")["reasons"], "ambiguous_base")
            self.assertEqual(r("no_init")["reasons"], "no_initiator")
            self.assertEqual(r("internal_stop")["reasons"], "internal_stop")
            self.assertEqual(r("no_stop")["reasons"], "no_stop")
            self.assertEqual(set(r("frame")["reasons"].split(",")), {"frame_length", "seq_frame_length", "no_stop"})
            self.assertEqual(r("phase")["reasons"], "phase_inconsistent")
            self.assertIn("overlapping_rows", r("overlap")["reasons"])
            self.assertNotIn(ids["pseudo"], rows)
            self.assertNotIn(ids["nc"], rows)
            self.assertNotIn(ids["iso_short"], rows)
            self.assertEqual(r("iso_long")["status"], "admitted")
            self.assertEqual(r("ov1")["reasons"], "topology")
            self.assertEqual(r("ov2")["reasons"], "topology")
            self.assertEqual(summary["topology_components"], 1)
            self.assertEqual(summary["benchmark_filtered"], {"biotype=lncRNA": 1, "pseudogene": 1})
            self.assertEqual(summary["nonrepresentative_omitted"], 1)
            self.assertEqual(summary["admitted"], 3)
            self.assertEqual(summary["masked_union"], len(rows) - 3)
            # every masked span is the full CDS span; the union counts bases once
            self.assertGreaterEqual(summary["masked_oriented_bases"],
                                    max(rw["cds_end"] - rw["cds_start"] + 1 for rw in rows.values() if rw["status"] == "masked"))
            # auxiliary sites: partial and sequence-incompatible ends are unknown, not negative
            self.assertGreaterEqual(summary["auxiliary_sites_unknown"]["start"], 3)   # edge partials + no_initiator
            self.assertGreaterEqual(summary["auxiliary_sites_unknown"]["stop"], 3)    # interior partial + no_stop + frame

    def test_prefix_marginalization_never_trims_exons(self):
        # a 5'-partial chain with phase 1 on the first row: two observed
        # bases of the partial codon, then in frame; the second exon keeps
        # every base (the proposal forbids discarding leading bases per exon)
        g = Genome()
        tid = g.gene("chrP", "+", 0, ["AGAAT", "AA"], [25], phases=[1, 2], attrs=";start_range=.,1")
        rows, summary = self.audit(g)
        self.assertEqual(rows[tid]["status"], "admitted")
        self.assertEqual(rows[tid]["missing_prefix"], 2)

    def test_partial_end_must_touch_its_own_edge(self):
        # engels-0039 finding 1: a declared partial end is admissible only at
        # the edge of that end in transcriptional orientation; the opposite
        # end touching the boundary does not count. One chain per genome so
        # the fixtures cannot overlap or overwrite each other's bases.
        for strand, seqid in (("+", "chrP"), ("-", "chrB")):
            L = len(Genome().seqs[seqid])
            cases = {
                # 5'-partial declared, but the chain sits at the 3' edge
                "wrong5": (L - 8, ["AAGCCTAA"], [2], five_partial(strand), "masked", "partial_5,partial_away_from_edge"),
                # 3'-partial declared, but the chain sits at the 5' edge
                "wrong3": (0, ["ATGAAA"], [0], three_partial(strand), "masked", "partial_3,partial_away_from_edge"),
                # the same declarations at their own edges are admitted
                "right5": (0, ["AAGCCTAA"], [2], five_partial(strand), "admitted", "partial_5"),
                "right3": (L - 6, ["ATGAAA"], [0], three_partial(strand), "admitted", "partial_3"),
                # both ends declared on a chain touching only one edge: masked
                "both_one": (0, ["AAGCCAAA"], [2], five_partial(strand) + three_partial(strand),
                             "masked", "partial_5,partial_3,partial_away_from_edge"),
                # both ends declared on a chain spanning the whole sequence: admitted
                "both": (0, ["AAGCC" + "C" * (L - 8) + "AAA"], [2], five_partial(strand) + three_partial(strand),
                         "admitted", "partial_5,partial_3"),
            }
            for name, (pos, pieces, phases, attrs, status, reasons) in cases.items():
                g = Genome()
                tid = g.gene(seqid, strand, pos, pieces, [], phases=phases, attrs=attrs)
                rows, _ = self.audit(g)
                self.assertEqual(rows[tid]["status"], status, (strand, name))
                self.assertEqual(rows[tid]["reasons"], reasons, (strand, name))

    def test_transcript_level_transl_except_masks(self):
        # engels-0039 finding 3: transl_except on the mRNA row alone
        for strand, seqid in (("+", "chrP"), ("-", "chrB")):
            g = Genome()
            tid = g.gene(seqid, strand, 100, ["ATGAAATAA"], [], tx_attrs=";transl_except=(pos:1..3,aa:Met)")
            rows, _ = self.audit(g)
            self.assertEqual(rows[tid]["status"], "masked", strand)
            self.assertEqual(rows[tid]["reasons"], "transl_except", strand)

    def test_auxiliary_sites_need_a_reliable_observation(self):
        # engels-0039 finding 2: alternative transcripts are sequence-audited
        # too; a start observed only by an initiator-less alternative is unknown
        for strand, seqid in (("+", "chrP"), ("-", "chrB")):
            g = Genome()
            g.gene(seqid, strand, 100, ["CTGAAATAA"], [], gene_id="shared")
            g.gene(seqid, strand, 400, ["ATGAAAAAATAA"], [], gene_id="shared")
            _, summary = self.audit(g)
            self.assertEqual(summary["auxiliary_sites_annotated"]["start"], 2, strand)
            self.assertEqual(summary["auxiliary_sites_retained"]["start"], 1, strand)
            self.assertEqual(summary["auxiliary_sites_unknown"]["start"], 1, strand)
            self.assertEqual(summary["auxiliary_sites_retained"]["stop"], 2, strand)
            # the same site observed reliably by another transcript stays positive
            g = Genome()
            g.gene(seqid, strand, 100, ["CTGAAATAA"], [], gene_id="shared")
            g.gene(seqid, strand, 100, ["ATGAAATAA"], [], gene_id="shared")
            _, summary = self.audit(g)
            self.assertEqual(summary["auxiliary_sites_annotated"]["start"], 1, strand)
            self.assertEqual(summary["auxiliary_sites_retained"]["start"], 1, strand)
            self.assertEqual(summary["auxiliary_sites_unknown"]["start"], 0, strand)

    def test_interior_partial_makes_its_junction_unknown(self):
        # engels-0039 finding 2, second case: a range declaration at the first
        # row's 3' boundary masks the chain and makes that donor unknown; the
        # acceptor across the intron keeps its observation
        for strand, seqid in (("+", "chrP"), ("-", "chrB")):
            g = Genome()
            tid = g.gene(seqid, strand, 100, ["ATGAAA", "TAA"], [25], attrs=three_partial(strand))
            rows, summary = self.audit(g)
            self.assertEqual(rows[tid]["status"], "masked", strand)
            self.assertIn("partial_unlocated", rows[tid]["reasons"].split(","), strand)
            self.assertNotIn("partial_3", rows[tid]["reasons"].split(","), strand)
            self.assertEqual(summary["auxiliary_sites_annotated"]["donor"], 1, strand)
            self.assertEqual(summary["auxiliary_sites_retained"]["donor"], 0, strand)
            self.assertEqual(summary["auxiliary_sites_unknown"]["donor"], 1, strand)
            self.assertEqual(summary["auxiliary_sites_retained"]["acceptor"], 1, strand)
            self.assertEqual(summary["auxiliary_sites_retained"]["start"], 1, strand)
            self.assertEqual(summary["auxiliary_sites_retained"]["stop"], 1, strand)
            # the mirror: a range at the second row's 5' boundary hides the acceptor
            g = Genome()
            tid = g.gene(seqid, strand, 100, ["ATGAAA", "TAA"], [25], row_attrs=["", five_partial(strand)])
            rows, summary = self.audit(g)
            self.assertIn("partial_unlocated", rows[tid]["reasons"].split(","), strand)
            self.assertEqual(summary["auxiliary_sites_retained"]["donor"], 1, strand)
            self.assertEqual(summary["auxiliary_sites_unknown"]["acceptor"], 1, strand)

    def test_repeated_transcript_ids_keep_their_own_audit(self):
        # stalin-0046: a transcript id reused on another sequence (or strand)
        # must never overwrite the other record's metadata, sequence result or
        # topology membership, whatever the FASTA order
        for strand in ("+", "-"):
            for order in (("chrP", "chrB"), ("chrB", "chrP")):
                g = Genome()
                g.gene("chrP", strand, 100, ["CTGAAATAA"], [], tid="shared")
                g.gene("chrB", strand, 100, ["ATGAAATAA"], [], tid="shared")
                rows, _ = self.audit_keyed(g, order)
                self.assertEqual(rows[("chrP", strand, "shared")]["reasons"], "no_initiator", (strand, order))
                self.assertEqual(rows[("chrB", strand, "shared")]["status"], "admitted", (strand, order))
            g = Genome()
            g.gene("chrP", strand, 100, ["ATGAAATAA"], [], tid="shared", attrs=";exception=ribosomal slippage")
            g.gene("chrB", strand, 100, ["ATGAAATAA"], [], tid="shared")
            rows, _ = self.audit_keyed(g)
            self.assertEqual(rows[("chrP", strand, "shared")]["reasons"], "exception", strand)
            self.assertEqual(rows[("chrB", strand, "shared")]["status"], "admitted", strand)
            # topology: only the overlapping chrP pair is masked, not chrB's copy
            g = Genome()
            g.gene("chrP", strand, 100, ["ATGAAA", "TAA"], [60], tid="shared")
            g.gene("chrB", strand, 100, ["ATGAAATAA"], [], tid="shared")
            g.gene("chrP", strand, 120, ["ATGCCCTAA"], [])
            rows, summary = self.audit_keyed(g)
            self.assertEqual(rows[("chrP", strand, "shared")]["reasons"], "topology", strand)
            self.assertEqual(rows[("chrB", strand, "shared")]["reasons"], "-", strand)
            self.assertEqual(summary["topology_masked"], 2, strand)
            self.assertEqual(sum(1 for r in rows.values() if "topology" in r["reasons"]), 2, strand)
        # the same id on both strands of one sequence: independent records
        g = Genome()
        g.gene("chrP", "+", 100, ["CTGAAATAA"], [], tid="shared")
        g.gene("chrP", "-", 100, ["ATGAAATAA"], [], tid="shared")
        rows, _ = self.audit_keyed(g)
        self.assertEqual(rows[("chrP", "+", "shared")]["reasons"], "no_initiator")
        self.assertEqual(rows[("chrP", "-", "shared")]["status"], "admitted")

    def test_table_conflict_and_alt_initiator(self):
        g = Genome()
        t1 = g.gene("chrP", "+", 100, ["ATGAAATAA"], [], attrs=";transl_table=6")
        rows, _ = self.audit(g)
        self.assertEqual(rows[t1]["reasons"], "table_conflict")


if __name__ == "__main__":
    unittest.main()
