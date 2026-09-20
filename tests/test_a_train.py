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

    def test_source_rejects_unknown_key(self):
        # A ``dev_seqid`` typo must not be silently dropped (stalin-0076).
        with self.assertRaises(ValueError):
            T.SpeciesSource.from_dict(self._src(dev_seqid=["chrDev"]))

    def test_source_rejects_string_dev_seqids(self):
        # ``list("chrDev")`` would reserve six one-letter ids; reject a string.
        with self.assertRaises(ValueError):
            T.SpeciesSource.from_dict(self._src(dev_seqids="chrDev"))

    def test_source_accepts_list_dev_seqids(self):
        s = T.SpeciesSource.from_dict(self._src(dev_seqids=["chrDev", "chr2"]))
        self.assertEqual(s.dev_seqids, ["chrDev", "chr2"])

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
    def _cfg(self, d, **kw):
        summary = os.path.join(d, "sp.summary.json")
        with open(summary, "w", encoding="utf-8") as fh:
            json.dump({"table": 6, "m": 30, "gff_md5": "aa", "fasta_md5": "bb"}, fh)
        kw.setdefault("out_dir", d)
        return T.TrainConfig(
            sources=[T.SpeciesSource("sp", summary, "g.gff", "f.fna",
                                     dev_seqids=["chrDev"])], **kw)

    def test_prefit_manifest_declares_plan_and_provenance(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = self._cfg(d, seed=3, steps=10, batch_size=4, eval_every=5)
            man = T.build_manifest(cfg, torch_version="x", cuda=None)
        self.assertEqual(man["task"], "T-human-014")
        self.assertEqual(man["scope"], "encoder-only-fixed-grammar")
        self.assertEqual(man["param_count"], SECTION_35_PARAM_COUNT)
        self.assertEqual(man["hyperparams"]["seed"], 3)
        self.assertEqual(man["hyperparams"]["eval_every"], 5)
        self.assertEqual(man["sampling_plan"]["planned_draws"], 40)  # steps*batch
        self.assertEqual(man["sources"][0]["name"], "sp")
        self.assertEqual(man["sources"][0]["dev_seqids"], ["chrDev"])
        self.assertEqual(man["sources"][0]["gff_md5"], "aa")
        self.assertEqual(man["sources"][0]["table"], 6)
        self.assertNotIn("actual", man)
        # The executed source is pinned independently of HEAD (engels-0080 P2).
        src = man["source"]
        self.assertEqual(len(src["source_sha256"]), 64)
        self.assertGreater(src["source_files"], 0)
        self.assertIn(src["source_dirty"], (True, False, None))

    def test_dirty_paths_keep_porcelain_status_columns(self):
        # engels-0081 P3: " M" (unstaged) has a blank first column; stripping
        # the line before slicing dropped the path's first character.
        from unittest.mock import patch
        for status, want in [
            (" M model/a/train.py\n", ["model/a/train.py"]),
            ("M  model/a/train.py\n", ["model/a/train.py"]),
            ("?? model/a/new.py\n", ["model/a/new.py"]),
            (" M model/a/train.py\nA  model/grammar/x.py\n",
             ["model/a/train.py", "model/grammar/x.py"]),
            ("", []),
        ]:
            with patch.object(T.subprocess, "check_output",
                              side_effect=[os.getcwd() + "\n", status]):
                got = T._source_provenance()
            self.assertEqual(got["source_dirty_files"], want, status)
            self.assertEqual(got["source_dirty"], bool(want), status)

    def test_source_digest_tracks_content(self):
        a = T._source_provenance()
        b = T._source_provenance()
        self.assertEqual(a["source_sha256"], b["source_sha256"])
        self.assertEqual(a["source_dirs"], ["model/a", "model/grammar"])

    def test_dev_seqids_and_eval_every_break_collision(self):
        # engels-0073: differing reserved chromosomes / eval schedule must yield
        # distinct manifests even at identical result counts.
        with tempfile.TemporaryDirectory() as d:
            summary = os.path.join(d, "sp.summary.json")
            with open(summary, "w", encoding="utf-8") as fh:
                json.dump({"table": 1, "m": 20, "gff_md5": "a", "fasta_md5": "b"}, fh)
            base = dict(out_dir=d, steps=10)
            a = T.TrainConfig(sources=[T.SpeciesSource(
                "sp", summary, "g", "f", dev_seqids=["chrA"])], eval_every=10, **base)
            b = T.TrainConfig(sources=[T.SpeciesSource(
                "sp", summary, "g", "f", dev_seqids=["chrB"])], eval_every=100, **base)
            ma = T.build_manifest(a, torch_version="x", cuda=None)
            mb = T.build_manifest(b, torch_version="x", cuda=None)
        self.assertNotEqual(ma, mb)

    def test_record_actual_appends_executed_work(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = self._cfg(d, seed=3, steps=10, batch_size=4)
            man = T.build_manifest(cfg, torch_version="x", cuda=None)
            man = T.record_actual(
                man, attempted_draws=40, accepted_windows=37, sampled_bases=1234,
                train_windows=12, dev_windows=3, best_dev_nll=0.5)
        self.assertEqual(man["actual"]["attempted_draws"], 40)
        self.assertEqual(man["actual"]["accepted_windows"], 37)
        self.assertEqual(man["actual"]["sampled_bases"], 1234)
        self.assertEqual(man["actual"]["train_windows"], 12)


class TestReservations(unittest.TestCase):
    @dataclass
    class FakeStats:
        windows_by_seqid: dict

    def _cfg(self, dev):
        return T.TrainConfig(
            sources=[T.SpeciesSource("sp", "s.json", "g", "f", dev_seqids=dev)],
            out_dir="/tmp/o")

    def test_declared_seqid_present_passes(self):
        stats = {"sp": self.FakeStats({"chrDev": 3, "chr1": 5})}
        T.validate_dev_reservations(self._cfg(["chrDev"]), stats)  # no raise

    def test_missing_declared_seqid_rejected(self):
        stats = {"sp": self.FakeStats({"chr1": 5})}
        with self.assertRaises(ValueError):
            T.validate_dev_reservations(self._cfg(["chrDev"]), stats)

    def test_no_reservation_passes(self):
        stats = {"sp": self.FakeStats({"chr1": 5})}
        T.validate_dev_reservations(self._cfg([]), stats)  # no raise


if __name__ == "__main__":
    unittest.main()
