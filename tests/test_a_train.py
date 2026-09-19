"""Torch-free unit tests for the candidate-A training entry point.

The tensor path (encoder forward, chain loss, autograd, checkpoint select) runs
only on a torch host and is exercised on gagarin; here we pin the pure-Python
scaffolding the fitting run depends on: config validation, window padding to the
pooling stride, the by-sequence train/dev split (the leakage rule), and the run
manifest. These need no torch, so they run on any interpreter.
"""
import json
import os
import tempfile
import unittest
from dataclasses import dataclass
from typing import Tuple

from model.a.inventory import POOL_STRIDE, SECTION_35_PARAM_COUNT
from model.a import train as T


@dataclass
class FakeWindow:
    key: Tuple[str, str, str]
    n: int


class TestConfig(unittest.TestCase):
    def _src(self, **kw):
        d = {"name": "sp", "summary": "s.json", "gff": "g.gff", "fasta": "f.fna"}
        d.update(kw)
        return d

    def test_source_requires_paths(self):
        with self.assertRaises(ValueError):
            T.SpeciesSource.from_dict({"name": "sp"})

    def test_config_parses_and_defaults(self):
        c = T.TrainConfig.from_dict({"sources": [self._src()], "out_dir": "/tmp/o"})
        self.assertEqual(len(c.sources), 1)
        self.assertEqual(c.sources[0].name, "sp")
        self.assertEqual(c.seed, 0)
        self.assertEqual(c.steps, 1000)
        self.assertIsNone(c.max_window)

    def test_config_overrides(self):
        c = T.TrainConfig.from_dict({
            "sources": [self._src(dev_seqids=["chrDev"])],
            "out_dir": "/tmp/o", "seed": 7, "steps": 5, "batch_size": 2,
            "lr": 1e-3, "max_window": 49152, "device": "cuda"})
        self.assertEqual(c.seed, 7)
        self.assertEqual(c.max_window, 49152)
        self.assertEqual(c.device, "cuda")
        self.assertEqual(c.sources[0].dev_seqids, ["chrDev"])

    def test_config_rejects_empty_sources(self):
        with self.assertRaises(ValueError):
            T.TrainConfig.from_dict({"sources": [], "out_dir": "/tmp/o"})

    def test_config_rejects_missing_out_dir(self):
        with self.assertRaises(ValueError):
            T.TrainConfig.from_dict({"sources": [self._src()]})

    def test_config_rejects_duplicate_species(self):
        with self.assertRaises(ValueError):
            T.TrainConfig.from_dict(
                {"sources": [self._src(), self._src()], "out_dir": "/tmp/o"})

    def test_config_rejects_unknown_keys(self):
        with self.assertRaises(ValueError):
            T.TrainConfig.from_dict(
                {"sources": [self._src()], "out_dir": "/tmp/o", "lrate": 1e-3})


class TestWindowing(unittest.TestCase):
    def test_pad_length_is_next_stride_multiple(self):
        self.assertEqual(T.pad_length(1), POOL_STRIDE)
        self.assertEqual(T.pad_length(POOL_STRIDE), POOL_STRIDE)
        self.assertEqual(T.pad_length(POOL_STRIDE + 1), 2 * POOL_STRIDE)

    def test_pad_length_rejects_nonpositive(self):
        with self.assertRaises(ValueError):
            T.pad_length(0)

    def test_pad_window_marks_padding_unavailable(self):
        seq = "ACGT"
        padded, avail = T.pad_window(seq)
        self.assertEqual(len(padded), T.pad_length(len(seq)))
        self.assertEqual(len(avail), len(padded))
        self.assertTrue(all(avail[:4]))
        self.assertFalse(any(avail[4:]))
        self.assertTrue(padded.startswith(seq))
        self.assertEqual(set(padded[4:]), {"N"})

    def test_pad_window_exact_multiple_unchanged(self):
        seq = "A" * POOL_STRIDE
        padded, avail = T.pad_window(seq)
        self.assertEqual(padded, seq)
        self.assertTrue(all(avail))


class TestSplit(unittest.TestCase):
    def test_split_by_seqid(self):
        ex = [
            FakeWindow(("chr1", "+", "t1"), 100),
            FakeWindow(("chrDev", "+", "t2"), 200),
            FakeWindow(("chr2", "-", "t3"), 300),
            FakeWindow(("chrDev", "-", "t4"), 400),
        ]
        train, dev = T.split_windows(ex, ["chrDev"])
        self.assertEqual({w.key[0] for w in train}, {"chr1", "chr2"})
        self.assertEqual({w.key[0] for w in dev}, {"chrDev"})
        self.assertEqual(len(dev), 2)

    def test_empty_dev_when_none_declared(self):
        ex = [FakeWindow(("chr1", "+", "t1"), 100)]
        train, dev = T.split_windows(ex, [])
        self.assertEqual(len(train), 1)
        self.assertEqual(dev, [])


class TestManifest(unittest.TestCase):
    def test_manifest_reads_digests_and_records_budget(self):
        with tempfile.TemporaryDirectory() as d:
            summary = os.path.join(d, "sp.summary.json")
            with open(summary, "w", encoding="utf-8") as fh:
                json.dump({"table": 6, "m": 30,
                           "gff_md5": "aa", "fasta_md5": "bb"}, fh)
            cfg = T.TrainConfig(
                sources=[T.SpeciesSource("sp", summary, "g.gff", "f.fna",
                                         dev_seqids=["chrDev"])],
                out_dir=d, seed=3, steps=10, batch_size=4)
            man = T.build_manifest(
                cfg, train_windows=12, dev_windows=3, sampled_bases=1234,
                best_dev_nll=0.5, torch_version="x", cuda=None)
        self.assertEqual(man["task"], "T-human-014")
        self.assertEqual(man["seed"], 3)
        self.assertEqual(man["param_count"], SECTION_35_PARAM_COUNT)
        self.assertEqual(man["sampled_bases"], 1234)
        self.assertEqual(man["train_windows"], 12)
        self.assertEqual(man["sources"][0]["gff_md5"], "aa")
        self.assertEqual(man["sources"][0]["table"], 6)
        self.assertEqual(man["sources"][0]["m"], 30)


if __name__ == "__main__":
    unittest.main()
