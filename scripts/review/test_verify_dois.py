#!/usr/bin/env python3
"""Regression tests for the DOI-verification TSV serialization.

The bug these cover: Crossref returns titles with embedded newlines (the
registered Helixer title breaks after "Gene"), and writing them verbatim
split one logical record across several physical lines, so `csv.reader`
saw 123 records of width {1, 3, 4, 6} for 111 bibliography entries and a
`reviews` value such as `marx` could be read as a `status`.

Standard library only:

    python3 scripts/review/test_verify_dois.py
"""

from __future__ import annotations

import csv
import io
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_dois import COLUMNS, check_serialization, flatten, serialize  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]

MULTILINE = (
    "holst2023helixer",
    "10.1093/nar/gkae1044",
    "ok",
    "Helixer—de novo Prediction of Primary Gene\nModels Combining Deep "
    "Learning\tand a Hidden Markov Model",
    "0.83",
    "engels, lenin, marx",
)


class TestSerialization(unittest.TestCase):
    def test_multiline_title_stays_one_record(self):
        text = serialize([MULTILINE])
        records = list(csv.reader(io.StringIO(text), delimiter="\t"))
        self.assertEqual(len(records), 2)
        self.assertEqual(len(records[1]), len(COLUMNS))
        self.assertEqual(records[1][0], "holst2023helixer")
        self.assertEqual(records[1][5], "engels, lenin, marx")
        self.assertNotIn("\n", records[1][3])

    def test_flatten_collapses_all_whitespace(self):
        self.assertEqual(flatten(" a\n\tb  c "), "a b c")

    def test_wrong_width_row_is_rejected(self):
        with self.assertRaises(ValueError):
            serialize([("key", "doi", "ok")])

    def test_check_rejects_a_split_record(self):
        broken = "\t".join(COLUMNS) + "\nk\td\tok\ttwo\nlines\t0.9\tlenin\n"
        with self.assertRaises(ValueError):
            check_serialization(broken, [MULTILINE])


class TestCommittedFile(unittest.TestCase):
    """The file in the tree must satisfy the same contract."""

    def setUp(self):
        path = ROOT / "docs" / "review" / "doi-verification.tsv"
        if not path.exists():
            self.skipTest(f"{path} not present")
        self.records = list(csv.reader(io.StringIO(path.read_text("utf-8")), delimiter="\t"))

    def test_header_and_widths(self):
        self.assertEqual(tuple(self.records[0]), COLUMNS)
        wrong = [i + 1 for i, r in enumerate(self.records) if len(r) != len(COLUMNS)]
        self.assertEqual(wrong, [], f"records with wrong width at lines {wrong}")

    def test_one_record_per_bibliography_entry(self):
        sys.path.insert(0, str(ROOT / "scripts" / "review"))
        from merge_refs import parse_bib  # noqa: E402

        bib = parse_bib((ROOT / "docs" / "refs" / "refs.bib").read_text("utf-8"))
        self.assertEqual([r[0] for r in self.records[1:]], [e["key"] for e in bib])


if __name__ == "__main__":
    unittest.main()
