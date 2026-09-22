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
        # chain-only scope by default: no gene-free windows
        self.assertEqual(c.background_windows, 0)
        self.assertEqual(c.background_length, 4096)

    def test_config_background_windows(self):
        c = T.TrainConfig.from_dict({"sources": [self._src()], "out_dir": "/tmp/o",
                                     "background_windows": 200, "background_length": 8192})
        self.assertEqual((c.background_windows, c.background_length), (200, 8192))
        with self.assertRaises(ValueError):
            T.TrainConfig.from_dict({"sources": [self._src()], "out_dir": "/tmp/o",
                                     "background_windows": -1})
        with self.assertRaises(ValueError):
            T.TrainConfig.from_dict({"sources": [self._src()], "out_dir": "/tmp/o",
                                     "background_length": 0})

    def test_config_context(self):
        c = T.TrainConfig.from_dict({"sources": [self._src()], "out_dir": "/tmp/o"})
        self.assertEqual(c.context, 0)
        c = T.TrainConfig.from_dict({"sources": [self._src()], "out_dir": "/tmp/o",
                                     "context": 512})
        self.assertEqual(c.context, 512)
        with self.assertRaises(ValueError):
            T.TrainConfig.from_dict({"sources": [self._src()], "out_dir": "/tmp/o",
                                     "context": -1})

    def test_config_background_length_bounded_by_max_window(self):
        # engels-0096: max_window bounds background windows too
        with self.assertRaises(ValueError):
            T.TrainConfig.from_dict({"sources": [self._src()], "out_dir": "/tmp/o",
                                     "max_window": 100, "background_windows": 1,
                                     "background_length": 1000})
        # equal length, a null maximum, or disabled background all pass
        T.TrainConfig.from_dict({"sources": [self._src()], "out_dir": "/tmp/o",
                                 "max_window": 1000, "background_windows": 1,
                                 "background_length": 1000})
        T.TrainConfig.from_dict({"sources": [self._src()], "out_dir": "/tmp/o",
                                 "background_windows": 1, "background_length": 1000})
        c = T.TrainConfig.from_dict({"sources": [self._src()], "out_dir": "/tmp/o",
                                     "max_window": 100, "background_length": 1000})
        # a directly constructed config is caught by the loader's own check
        c.background_windows = 1
        with self.assertRaises(ValueError):
            c.check_window_bound()

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
        self.assertEqual(man["scope"], "encoder-and-pooled-decoder")
        self.assertEqual(man["param_count"], SECTION_35_PARAM_COUNT)
        self.assertEqual(man["hyperparams"]["seed"], 3)
        self.assertEqual(man["hyperparams"]["eval_every"], 5)
        self.assertEqual(man["sampling_plan"]["planned_draws"], 40)  # steps*batch
        self.assertEqual(man["hyperparams"]["context"], 0)
        self.assertEqual(man["sampling_plan"]["context_per_side"], 0)
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
        # Optional cost fields default to empty so older callers still work.
        self.assertIsNone(man["actual"]["best_step"])
        self.assertEqual(man["actual"]["history"], [])
        self.assertEqual(man["actual"]["timing"], {})
        self.assertEqual(man["actual"]["composition"], {})

    def test_record_actual_keeps_composition(self):
        comp = {"cds_bases": 100, "intron_bases": 50, "intergenic_bases": 30,
                "background_draws": 2, "context_bases": 20}
        with tempfile.TemporaryDirectory() as d:
            cfg = self._cfg(d, seed=3, steps=10, batch_size=4)
            man = T.build_manifest(cfg, torch_version="x", cuda=None)
            man = T.record_actual(
                man, attempted_draws=4, accepted_windows=4, sampled_bases=180,
                train_windows=12, dev_windows=3, best_dev_nll=0.5, composition=comp)
        self.assertEqual(man["actual"]["composition"], comp)
        json.dumps(man)

    def test_record_actual_keeps_best_step_history_and_timing(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = self._cfg(d, seed=3, steps=10, batch_size=4)
            man = T.build_manifest(cfg, torch_version="x", cuda=None)
            hist = [{"step": 5, "train_nll": 2.0, "dev_nll": 3.0, "elapsed_s": 1.0},
                    {"step": 10, "train_nll": 1.0, "dev_nll": 2.5, "elapsed_s": 2.0}]
            man = T.record_actual(
                man, attempted_draws=40, accepted_windows=40, sampled_bases=1,
                train_windows=12, dev_windows=3, best_dev_nll=2.5, best_step=10,
                history=hist, timing={"fit_wall_s": 2.0, "fit_cpu_s": 1.9})
        self.assertEqual(man["actual"]["best_step"], 10)
        self.assertEqual(man["actual"]["history"], hist)
        self.assertEqual(man["actual"]["timing"]["fit_cpu_s"], 1.9)
        json.dumps(man)  # the manifest stays serialisable


class TestDevSubsample(unittest.TestCase):
    def test_none_or_large_limit_keeps_all_in_order(self):
        dev = ["a", "b", "c"]
        self.assertEqual(T.subsample_dev(dev, None, 0), dev)
        self.assertEqual(T.subsample_dev(dev, 3, 0), dev)
        self.assertEqual(T.subsample_dev(dev, 99, 0), dev)

    def test_limit_draws_without_replacement_in_load_order_and_is_seeded(self):
        dev = list(range(100))
        a = T.subsample_dev(dev, 10, 7)
        self.assertEqual(len(a), 10)
        self.assertEqual(len(set(a)), 10)
        self.assertEqual(a, sorted(a))
        self.assertEqual(a, T.subsample_dev(dev, 10, 7))
        self.assertNotEqual(a, T.subsample_dev(dev, 10, 8))

    def test_config_key_is_validated_and_recorded(self):
        with tempfile.TemporaryDirectory() as d:
            summary = os.path.join(d, "sp.summary.json")
            with open(summary, "w", encoding="utf-8") as fh:
                json.dump({"table": 6, "m": 30, "gff_md5": "aa", "fasta_md5": "bb"}, fh)
            base = {"sources": [{"name": "sp", "summary": summary, "gff": "g", "fasta": "f"}],
                    "out_dir": d}
            self.assertIsNone(T.TrainConfig.from_dict(base).dev_windows_max)
            cfg = T.TrainConfig.from_dict(dict(base, dev_windows_max=256))
            self.assertEqual(cfg.dev_windows_max, 256)
            man = T.build_manifest(cfg, torch_version="x", cuda=None)
            self.assertEqual(man["hyperparams"]["dev_windows_max"], 256)
            self.assertEqual(man["sampling_plan"]["dev_windows_max"], 256)
            with self.assertRaises(ValueError):
                T.TrainConfig.from_dict(dict(base, dev_windows_max=0))
            man = T.record_actual(
                man, attempted_draws=1, accepted_windows=1, sampled_bases=1,
                train_windows=1, dev_windows=4893, best_dev_nll=None,
                dev_windows_evaluated=256)
            self.assertEqual(man["actual"]["dev_windows"], 4893)
            self.assertEqual(man["actual"]["dev_windows_evaluated"], 256)


class TestLearningRateSchedule(unittest.TestCase):
    """The schedule of ``lr_at`` (fit v3 runs "cosine"; v1/v2 ran "constant")."""

    def _cfg(self, **kw):
        kw.setdefault("steps", 10)
        return T.TrainConfig(
            sources=[T.SpeciesSource("sp", "s.json", "g", "f")],
            out_dir="/tmp/o", lr=1.0, **kw)

    def test_constant_holds_lr_at_every_step(self):
        cfg = self._cfg()
        self.assertEqual([T.lr_at(i, cfg) for i in range(10)], [1.0] * 10)

    def test_constant_ignores_warmup_and_floor(self):
        cfg = self._cfg(warmup_steps=4, lr_min_factor=0.5)
        self.assertEqual([T.lr_at(i, cfg) for i in range(10)], [1.0] * 10)

    def test_cosine_warmup_ramps_linearly_and_is_never_zero(self):
        cfg = self._cfg(lr_schedule="cosine", warmup_steps=4)
        self.assertEqual([T.lr_at(i, cfg) for i in range(4)],
                         [0.25, 0.5, 0.75, 1.0])

    def test_cosine_without_warmup_starts_at_lr(self):
        cfg = self._cfg(lr_schedule="cosine")
        self.assertEqual(T.lr_at(0, cfg), 1.0)

    def test_cosine_decays_monotonically_to_the_floor(self):
        cfg = self._cfg(lr_schedule="cosine", warmup_steps=2, lr_min_factor=0.05,
                        steps=100)
        tail = [T.lr_at(i, cfg) for i in range(2, 100)]
        self.assertEqual(tail, sorted(tail, reverse=True))
        self.assertAlmostEqual(tail[0], 1.0)
        self.assertAlmostEqual(tail[-1], 0.05)

    def test_cosine_floor_zero_ends_at_zero(self):
        cfg = self._cfg(lr_schedule="cosine", steps=50)
        self.assertAlmostEqual(T.lr_at(49, cfg), 0.0)

    def test_steps_past_the_end_stay_at_the_floor(self):
        cfg = self._cfg(lr_schedule="cosine", lr_min_factor=0.1)
        self.assertAlmostEqual(T.lr_at(50, cfg), 0.1)

    def test_config_validates_and_manifest_records_the_schedule(self):
        with tempfile.TemporaryDirectory() as d:
            summary = os.path.join(d, "sp.summary.json")
            with open(summary, "w", encoding="utf-8") as fh:
                json.dump({"table": 6, "m": 30, "gff_md5": "aa", "fasta_md5": "bb"}, fh)
            base = {"sources": [{"name": "sp", "summary": summary, "gff": "g", "fasta": "f"}],
                    "out_dir": d, "steps": 100}
            default = T.TrainConfig.from_dict(base)
            self.assertEqual(default.lr_schedule, "constant")
            self.assertEqual(default.warmup_steps, 0)
            self.assertEqual(default.lr_min_factor, 0.0)
            cfg = T.TrainConfig.from_dict(dict(base, lr_schedule="cosine",
                                               warmup_steps=10, lr_min_factor=0.05))
            man = T.build_manifest(cfg, torch_version="x", cuda=None)
            self.assertEqual(man["hyperparams"]["lr_schedule"], "cosine")
            self.assertEqual(man["hyperparams"]["warmup_steps"], 10)
            self.assertEqual(man["hyperparams"]["lr_min_factor"], 0.05)
            for bad in ({"lr_schedule": "linear"}, {"warmup_steps": -1},
                        {"warmup_steps": 100}, {"lr_min_factor": 1.5},
                        {"lr_min_factor": -0.1}):
                with self.assertRaises(ValueError):
                    T.TrainConfig.from_dict(dict(base, **bad))

    def test_cosine_needs_two_post_warmup_steps(self):
        """engels-0105/stalin-0110 P3: a one-step decay cannot both start at
        ``lr`` and end at the floor, so ``from_dict`` rejects it rather than
        silently returning the full rate at the final step. The rejected set
        is exactly ``warmup_steps == steps - 1`` (and ``steps=1``); constant
        schedules, which have no endpoint contract, are unaffected."""
        with tempfile.TemporaryDirectory() as d:
            summary = os.path.join(d, "sp.summary.json")
            with open(summary, "w", encoding="utf-8") as fh:
                json.dump({"table": 6, "m": 30, "gff_md5": "aa", "fasta_md5": "bb"}, fh)
            base = {"sources": [{"name": "sp", "summary": summary, "gff": "g", "fasta": "f"}],
                    "out_dir": d, "lr_schedule": "cosine"}
            for steps, warmup in ((1, 0), (10, 9), (3000, 2999)):
                with self.assertRaises(ValueError):
                    T.TrainConfig.from_dict(
                        dict(base, steps=steps, warmup_steps=warmup))
            for steps, warmup in ((2, 0), (10, 8), (3000, 150)):
                cfg = T.TrainConfig.from_dict(
                    dict(base, steps=steps, warmup_steps=warmup,
                         lr_min_factor=0.05, lr=3e-4))
                self.assertAlmostEqual(T.lr_at(steps - 1, cfg), 3e-4 * 0.05)
            # A constant schedule keeps every accepted (steps, warmup) pair.
            T.TrainConfig.from_dict(
                dict(base, lr_schedule="constant", steps=10, warmup_steps=9))

    def test_every_accepted_cosine_endpoint_reaches_the_floor(self):
        """stalin-0110's sweep, bounded: over ``steps`` 1..64 and every legal
        ``warmup_steps``, each config that ``from_dict`` accepts ends its
        final step at the floor."""
        with tempfile.TemporaryDirectory() as d:
            summary = os.path.join(d, "sp.summary.json")
            with open(summary, "w", encoding="utf-8") as fh:
                json.dump({"table": 6, "m": 30, "gff_md5": "aa", "fasta_md5": "bb"}, fh)
            base = {"sources": [{"name": "sp", "summary": summary, "gff": "g", "fasta": "f"}],
                    "out_dir": d, "lr_schedule": "cosine", "lr": 3e-4,
                    "lr_min_factor": 0.05}
            accepted = rejected = 0
            for steps in range(1, 65):
                for warmup in range(0, steps):
                    try:
                        cfg = T.TrainConfig.from_dict(
                            dict(base, steps=steps, warmup_steps=warmup))
                    except ValueError:
                        rejected += 1
                        self.assertEqual(warmup, steps - 1)
                        continue
                    accepted += 1
                    self.assertAlmostEqual(T.lr_at(steps - 1, cfg), 3e-4 * 0.05,
                                           delta=1e-15)
            self.assertEqual(rejected, 64)
            self.assertEqual(accepted, 2080 - 64)


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
