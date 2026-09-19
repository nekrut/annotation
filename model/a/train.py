"""Candidate A training and end-to-end measurement entry point (T-human-014).

This wires the three reviewed pieces of candidate A into one runnable program:

  * the compact DNA encoder (:mod:`model.a.encoder`, proposal section 3.5),
  * the structured training-window loader (:mod:`model.a.dataset`,
    section 3.6), which hard-verifies each species' pinned source digests and
    delegates admission to ``model.labels``, and
  * the differentiable chain loss (:mod:`model.a.torch_loss`), the
    reference-recurrence ``log Z - log Z_num`` whose gradient is the
    free-minus-numerator emission posterior.

Two subcommands:

  ``train``    fit the encoder + pooled-decoder scalars on the *train*
               sequences of the configured train species, selecting the
               checkpoint on the declared development sequences only, and
               write a run manifest (data digests, seed, commit, hardware,
               sampled bases) next to the checkpoint.

  ``measure``  time preprocessing, encoder and decoder/traceback separately
               over one sequence, with peak host and (if present) device
               memory, and emit a row in the ``docs/cost-baseline`` TSV
               convention for ``docs/design/a-pilot.md``.

Torch is imported lazily inside the functions that need a tensor, so the
config parsing, the train/dev split, the sampled-base accounting and the
manifest construction are all importable and unit-testable without torch on
the host. Fitting and measurement themselves run on gagarin under the Phase 4
budget (``relay/TASK.md``); nothing here trains a model on import.

The loss forward uses the *reference* recurrence (explicit mandatory-intron
states, Python-loop over boundaries), the same one the oracle and the
torch-parity tests pin. The fast vectorized delayed-entry kernel is the
follow-on; this entry point exists so the pilot can be measured end to end and
the tensor path verified on a real torch host, which no local run can do.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import resource
import subprocess
import sys
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .inventory import POOL_STRIDE, SECTION_35_PARAM_COUNT

# Ranges are half-open oriented intervals, the loader/loss convention.
Range = Tuple[int, int]


# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
@dataclass
class SpeciesSource:
    """One pinned train species and its train/dev sequence split.

    ``summary``/``gff``/``fasta`` are the paths the loader hard-verifies against
    the committed digests. ``dev_seqids`` are the sequence ids reserved for
    checkpoint selection and development scoring; every other admitted sequence
    is a training sequence. The split is by *sequence* (chromosome/scaffold),
    the leakage rule Phase 4 requires: a dev sequence never contributes a
    gradient step.
    """

    name: str
    summary: str
    gff: str
    fasta: str
    dev_seqids: List[str] = field(default_factory=list)

    @classmethod
    def from_dict(cls, d: dict) -> "SpeciesSource":
        missing = {"name", "summary", "gff", "fasta"} - set(d)
        if missing:
            raise ValueError(f"species source missing keys: {sorted(missing)}")
        return cls(
            name=str(d["name"]),
            summary=str(d["summary"]),
            gff=str(d["gff"]),
            fasta=str(d["fasta"]),
            dev_seqids=list(d.get("dev_seqids", [])),
        )


@dataclass
class TrainConfig:
    """A fully declared fitting run. Everything a manifest must record."""

    sources: List[SpeciesSource]
    out_dir: str
    seed: int = 0
    steps: int = 1000
    batch_size: int = 8
    lr: float = 3e-4
    weight_decay: float = 0.0
    grad_clip: float = 1.0
    max_window: Optional[int] = None
    eval_every: int = 100
    device: str = "cpu"

    @classmethod
    def from_dict(cls, d: dict) -> "TrainConfig":
        if "sources" not in d or not d["sources"]:
            raise ValueError("config needs a non-empty 'sources' list")
        if "out_dir" not in d:
            raise ValueError("config needs 'out_dir'")
        sources = [SpeciesSource.from_dict(s) for s in d["sources"]]
        names = [s.name for s in sources]
        if len(set(names)) != len(names):
            raise ValueError("duplicate species names in config")
        known = {f.name for f in cls.__dataclass_fields__.values()}  # type: ignore[attr-defined]
        unknown = set(d) - known - {"sources"}
        if unknown:
            raise ValueError(f"unknown config keys: {sorted(unknown)}")
        kwargs = {k: d[k] for k in d if k in known and k != "sources"}
        return cls(sources=sources, **kwargs)

    @classmethod
    def from_json(cls, path: str) -> "TrainConfig":
        with open(path, "r", encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))


# --------------------------------------------------------------------------
# Windowing helpers (torch-free, so they are unit-testable without torch)
# --------------------------------------------------------------------------
def pad_length(n: int, stride: int = POOL_STRIDE) -> int:
    """Smallest multiple of ``stride`` that is >= ``n`` (the encoder needs
    ``L % POOL_STRIDE == 0``)."""
    if n <= 0:
        raise ValueError("window length must be positive")
    return ((n + stride - 1) // stride) * stride


def pad_window(window: str, stride: int = POOL_STRIDE) -> Tuple[str, List[bool]]:
    """Right-pad ``window`` with unavailable ``N`` bases to a stride multiple.

    Returns the padded string and the per-position availability mask the
    featurizer uses; padding invents no sequence (availability ``False``), so a
    padded position contributes nothing and never contaminates a real base's GC
    window. The real length is ``len(window)``; emissions are cropped back to it
    before the loss, whose ``_check_input`` enforces ``shape == (11, len(x))``.
    """
    n = len(window)
    target = pad_length(n, stride)
    available = [True] * n + [False] * (target - n)
    return window + ("N" * (target - n)), available


def split_windows(examples, dev_seqids) -> Tuple[list, list]:
    """Partition loaded windows into (train, dev) by source sequence id.

    ``WindowExample.key`` is ``(seqid, strand, transcript)``; a window whose
    seqid is in ``dev_seqids`` is a development window, used only for
    checkpoint selection, never for a gradient step.
    """
    dev_set = set(dev_seqids)
    train, dev = [], []
    for ex in examples:
        (dev if ex.key[0] in dev_set else train).append(ex)
    return train, dev


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True,
            stderr=subprocess.DEVNULL).strip()
    except Exception:  # not in a checkout, or git absent
        return "unknown"


def _source_digests(source: SpeciesSource) -> dict:
    """The pinned gff/fasta MD5s the loader verifies, read from the summary so
    the manifest records exactly what was fitted on."""
    with open(source.summary, "r", encoding="utf-8") as fh:
        meta = json.load(fh)
    return {
        "summary": source.summary,
        "table": int(meta.get("table", 1)),
        "m": int(meta.get("m", 20)),
        "gff_md5": meta.get("gff_md5"),
        "fasta_md5": meta.get("fasta_md5"),
    }


def build_manifest(config: TrainConfig, *, train_windows: int, dev_windows: int,
                   sampled_bases: int, best_dev_nll: Optional[float],
                   torch_version: str, cuda: Optional[str]) -> dict:
    """The run manifest: everything needed to reproduce and to audit the budget.

    Recorded in the task log and committed next to the checkpoint (the
    checkpoint itself is not committed; it exceeds nothing but lives in the
    gagarin scratch/return path)."""
    return {
        "task": "T-human-014",
        "commit": _git_commit(),
        "seed": config.seed,
        "param_count": SECTION_35_PARAM_COUNT,
        "hardware": platform.platform(),
        "python": platform.python_version(),
        "torch": torch_version,
        "cuda": cuda,
        "device": config.device,
        "steps": config.steps,
        "batch_size": config.batch_size,
        "lr": config.lr,
        "weight_decay": config.weight_decay,
        "grad_clip": config.grad_clip,
        "max_window": config.max_window,
        "train_windows": train_windows,
        "dev_windows": dev_windows,
        "sampled_bases": sampled_bases,
        "best_dev_nll": best_dev_nll,
        "sources": [_source_digests(s) for s in config.sources],
    }


# --------------------------------------------------------------------------
# Torch-dependent training / measurement (imported lazily)
# --------------------------------------------------------------------------
def _load_all_windows(config: TrainConfig):
    """Load every admitted clean complete-target window for the configured
    species, tagging each with its species genetic-code table."""
    from .dataset import iter_windows, LoaderStats

    examples = []
    stats_by_species: Dict[str, LoaderStats] = {}
    for src in config.sources:
        stats = LoaderStats()
        for ex in iter_windows(src.summary, src.gff, src.fasta,
                               max_window=config.max_window, stats=stats):
            examples.append(ex)
        stats_by_species[src.name] = stats
    return examples, stats_by_species


def _emissions_for(model, ex, device, dtype):
    """Featurize one window, run the encoder, and return the ``(11, n)`` torch
    emissions cropped to the real window length, ready for ``chain_nll``."""
    import torch

    from .features import encode_sequence

    padded, available = pad_window(ex.window)
    x = encode_sequence(padded, available=available).unsqueeze(0).to(device)  # [1,8,L]
    emissions = model.encoder(x)[0]                                           # [11,L]
    return emissions[:, : ex.n].to(dtype)


def _window_loss(model, ex, device, dtype):
    """Differentiable chain NLL for one window, or ``None`` if the numerator
    support admits no legal path (an unusable crop; section 3.6 says drop it)."""
    import torch

    from .torch_loss import chain_nll
    from model.grammar.codes import TABLES

    emissions = _emissions_for(model, ex, device, dtype)
    code = TABLES[ex.table]
    loss = chain_nll(ex.window, emissions, ex.cds_ranges, ex.intron_ranges, code=code)
    if not torch.isfinite(loss):
        return None
    return loss


def train(config: TrainConfig, log=print) -> dict:
    """Fit candidate A and return the run manifest.

    One gradient step accumulates ``batch_size`` per-window losses (windows vary
    in length, so they are summed rather than tensor-batched; the fast kernel
    will collate). Development NLL is evaluated every ``eval_every`` steps and
    the lowest-NLL checkpoint is kept. The manifest is written to
    ``out_dir/run_manifest.json`` and the checkpoint to ``out_dir/best.pt``.
    """
    import os

    import torch

    from .encoder import CandidateA

    torch.manual_seed(config.seed)
    device = torch.device(config.device)
    dtype = torch.float64  # matches the oracle/parity tests

    examples, stats = _load_all_windows(config)
    dev_ids = [i for s in config.sources for i in s.dev_seqids]
    train_ex, dev_ex = split_windows(examples, dev_ids)
    if not train_ex:
        raise ValueError("no training windows after the dev split")
    log(f"loaded {len(examples)} windows: {len(train_ex)} train, {len(dev_ex)} dev")

    model = CandidateA().to(device)
    assert model.num_parameters() == SECTION_35_PARAM_COUNT, model.num_parameters()
    opt = torch.optim.Adam(model.parameters(), lr=config.lr,
                           weight_decay=config.weight_decay)
    gen = torch.Generator().manual_seed(config.seed)

    os.makedirs(config.out_dir, exist_ok=True)
    best_dev = None
    sampled_bases = 0

    def evaluate() -> Optional[float]:
        if not dev_ex:
            return None
        model.eval()
        total, count = 0.0, 0
        with torch.no_grad():
            for ex in dev_ex:
                loss = _window_loss(model, ex, device, dtype)
                if loss is not None:
                    total += float(loss)
                    count += 1
        model.train()
        return total / count if count else None

    model.train()
    for step in range(config.steps):
        idx = torch.randint(len(train_ex), (config.batch_size,), generator=gen)
        opt.zero_grad()
        batch_loss, used = 0.0, 0
        for j in idx.tolist():
            ex = train_ex[j]
            loss = _window_loss(model, ex, device, dtype)
            if loss is None:
                continue
            (loss / config.batch_size).backward()
            batch_loss += float(loss)
            sampled_bases += ex.n
            used += 1
        if used:
            torch.nn.utils.clip_grad_norm_(model.parameters(), config.grad_clip)
            opt.step()
        if (step + 1) % config.eval_every == 0 or step + 1 == config.steps:
            dev_nll = evaluate()
            log(f"step {step + 1}: train_nll/window={batch_loss / max(used, 1):.4f} "
                f"dev_nll={dev_nll}")
            if dev_nll is not None and (best_dev is None or dev_nll < best_dev):
                best_dev = dev_nll
                torch.save(model.state_dict(), os.path.join(config.out_dir, "best.pt"))
    if best_dev is None:  # no dev set: keep the final model
        torch.save(model.state_dict(), os.path.join(config.out_dir, "best.pt"))

    manifest = build_manifest(
        config, train_windows=len(train_ex), dev_windows=len(dev_ex),
        sampled_bases=sampled_bases, best_dev_nll=best_dev,
        torch_version=torch.__version__,
        cuda=(torch.version.cuda if torch.cuda.is_available() else None))
    with open(os.path.join(config.out_dir, "run_manifest.json"), "w",
              encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
    return manifest


def _scores_from_emissions(emissions):
    """A CPU ``model.grammar.scores.Scores`` view of a ``(11, n)`` tensor, for
    the decoder/traceback stage of measurement (the reference decoder consumes
    Python-float channels). The tensor->list conversion is timed as part of the
    decode stage, not the encoder."""
    from model.grammar.scores import Scores

    e = emissions.detach().to("cpu").double().tolist()
    return Scores(n=emissions.shape[1], u=e[0], cds=[e[1], e[2], e[3]],
                  intron=[e[4], e[5], e[6]], start=e[7], stop=e[8],
                  donor=e[9], acceptor=e[10])


def measure(config: TrainConfig, species: str, seqid: Optional[str],
            checkpoint: Optional[str], log=print) -> dict:
    """Time preprocessing, encoder and decode/traceback over one sequence.

    Returns a dict of measured stage times, oriented bases, and per-Mb rates,
    plus peak host RSS and device memory. The CPU regime reports process CPU
    seconds per Mb; the GPU regime the wall seconds per Mb, matching
    ``docs/cost-baseline`` conventions. Windows on ``seqid`` (or every window if
    ``seqid`` is None) are processed; the decoder is the reference delayed-entry
    Viterbi over the emissions.
    """
    import torch

    from .encoder import CandidateA
    from model.grammar.codes import TABLES
    from model.grammar.delayed import DelayedEntryDecoder

    src = next((s for s in config.sources if s.name == species), None)
    if src is None:
        raise ValueError(f"species {species!r} not in config")
    device = torch.device(config.device)
    dtype = torch.float64

    model = CandidateA().to(device)
    if checkpoint:
        model.load_state_dict(torch.load(checkpoint, map_location=device))
    model.eval()

    examples, _ = _load_all_windows(
        TrainConfig(sources=[src], out_dir=config.out_dir,
                    max_window=config.max_window))
    if seqid is not None:
        examples = [e for e in examples if e.key[0] == seqid]
    if not examples:
        raise ValueError(f"no windows for {species} {seqid}")

    from .features import encode_sequence

    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)

    t_pre = t_enc = t_dec = 0.0
    bases = 0
    with torch.no_grad():
        for ex in examples:
            padded, available = pad_window(ex.window)
            c0 = time.process_time()
            feats = encode_sequence(padded, available=available).unsqueeze(0).to(device)
            t_pre += time.process_time() - c0

            c0 = time.process_time()
            emissions = model.encoder(feats)[0][:, : ex.n].to(dtype)
            if device.type == "cuda":
                torch.cuda.synchronize(device)
            t_enc += time.process_time() - c0

            c0 = time.process_time()
            decoder = DelayedEntryDecoder(TABLES[ex.table])
            sc = _scores_from_emissions(emissions)
            decoder.viterbi(ex.window, sc)
            t_dec += time.process_time() - c0
            bases += ex.n

    mb = bases / 1e6
    device_mem_gb = (torch.cuda.max_memory_allocated(device) / 2**30
                     if device.type == "cuda" else None)
    row = {
        "species": species, "seqid": seqid, "windows": len(examples),
        "oriented_bases": bases, "device": config.device,
        "preprocess_cpu_s": t_pre, "encoder_cpu_s": t_enc, "decode_cpu_s": t_dec,
        "cpu_s_per_mb": (t_pre + t_enc + t_dec) / mb if mb else None,
        "peak_host_rss_gb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20,
        "peak_device_mem_gb": device_mem_gb,
    }
    log(json.dumps(row, indent=2, sort_keys=True))
    return row


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    pt = sub.add_parser("train", help="fit candidate A and write a run manifest")
    pt.add_argument("--config", required=True, help="path to a TrainConfig JSON")

    pm = sub.add_parser("measure", help="time preprocessing/encoder/decode on one sequence")
    pm.add_argument("--config", required=True, help="path to a TrainConfig JSON")
    pm.add_argument("--species", required=True)
    pm.add_argument("--seqid", default=None, help="restrict to one sequence id")
    pm.add_argument("--checkpoint", default=None, help="a best.pt to load")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    config = TrainConfig.from_json(args.config)
    if args.cmd == "train":
        train(config)
    elif args.cmd == "measure":
        measure(config, args.species, args.seqid, args.checkpoint)
    return 0


if __name__ == "__main__":
    sys.exit(main())
