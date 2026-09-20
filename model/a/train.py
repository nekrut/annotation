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

  ``train``    fit the encoder against the fixed-grammar reference chain loss
               on the *train* sequences of the configured train species,
               selecting the checkpoint on the declared development sequences
               only, and write a run manifest (declared draw plan, data
               digests, seed, commit, hardware, actual sampled bases) next to
               the checkpoint. This is an encoder-only profiling fit: the
               pooled-decoder scalars have no gradient under the fixed grammar
               and are wired in a later increment (engels-0073).

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
from pathlib import Path
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
        known = {"name", "summary", "gff", "fasta", "dev_seqids"}
        unknown = set(d) - known
        if unknown:
            # A silently-ignored key (e.g. the ``dev_seqid`` typo) would drop a
            # declared reservation and train on the intended dev sequence
            # (stalin-0076); reject it up front.
            raise ValueError(f"unknown species source keys: {sorted(unknown)}")
        missing = {"name", "summary", "gff", "fasta"} - set(d)
        if missing:
            raise ValueError(f"species source missing keys: {sorted(missing)}")
        dev = d.get("dev_seqids", [])
        # A bare string is iterable per character; ``list("chrDev")`` would
        # reserve six one-letter sequences and leave the real chromosome in the
        # gradient pool. Require an explicit list/tuple of ids.
        if isinstance(dev, str) or not isinstance(dev, (list, tuple)):
            raise ValueError("dev_seqids must be a list of sequence ids")
        return cls(
            name=str(d["name"]),
            summary=str(d["summary"]),
            gff=str(d["gff"]),
            fasta=str(d["fasta"]),
            dev_seqids=[str(x) for x in dev],
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
    # "fast": the vectorized delayed-entry kernel (model.a.fast_loss);
    # "reference": the expanded-grammar torch forward (model.a.torch_loss),
    # 40-80x slower and ~10x the memory, kept for parity checks only.
    loss_kernel: str = "fast"

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
        if kwargs.get("loss_kernel", "fast") not in ("fast", "reference"):
            raise ValueError("loss_kernel must be 'fast' or 'reference'")
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


SOURCE_DIRS = ("model/a", "model/grammar")


def _source_provenance() -> dict:
    """Identify the *executed* source, not only ``HEAD``.

    ``commit`` alone cannot reproduce a run made from a dirty tree (the fast-
    kernel smoke run of 2026-09-20 was profiled before its source was committed;
    engels-0080). So the manifest also records whether ``git status`` reports
    changes under :data:`SOURCE_DIRS`, and a SHA-256 over the sorted contents
    of every ``*.py`` under them, which pins the code that ran whatever the
    tree's git state was.
    """
    root = None
    try:
        root = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], text=True,
            stderr=subprocess.DEVNULL).strip()
        dirty = subprocess.check_output(
            ["git", "status", "--porcelain", "--", *SOURCE_DIRS], text=True,
            stderr=subprocess.DEVNULL)
        # Columns 1-2 are status codes and may be blank (" M" = unstaged
        # modification); do not strip the line before slicing (engels-0081).
        dirty_files = sorted(line[3:] for line in dirty.splitlines()
                             if line.strip())
    except Exception:  # not in a checkout, or git absent
        dirty_files = None
    if root is None:
        root = str(Path(__file__).resolve().parents[2])
    digest = hashlib.sha256()
    files = []
    for d in SOURCE_DIRS:
        for f in sorted(Path(root, d).glob("*.py")):
            files.append(str(f.relative_to(root)))
            digest.update(files[-1].encode() + b"\0")
            digest.update(f.read_bytes() + b"\0")
    return {
        "source_dirs": list(SOURCE_DIRS),
        "source_files": len(files),
        "source_sha256": digest.hexdigest(),
        "source_dirty": None if dirty_files is None else bool(dirty_files),
        "source_dirty_files": dirty_files,
    }


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


def validate_dev_reservations(config: TrainConfig, stats_by_species) -> None:
    """Reject a declared development reservation that no admitted window can
    satisfy, *before* any gradient step.

    Each source's ``dev_seqids`` must name sequences that actually produced
    admitted windows for that species. A typo, a wrong id, or a reservation
    emptied by filtering would otherwise leave the intended development
    sequence in the gradient pool and let ``evaluate`` fall back to the final
    checkpoint (stalin-0076). ``stats_by_species[name].windows_by_seqid`` is the
    per-species admitted inventory the loader already tallies.
    """
    for s in config.sources:
        present = set(stats_by_species[s.name].windows_by_seqid)
        for seqid in s.dev_seqids:
            if seqid not in present:
                raise ValueError(
                    f"species {s.name!r} reserves dev seqid {seqid!r}, but no "
                    f"admitted window has that sequence id "
                    f"(present: {sorted(present)}); a declared reservation must "
                    f"match the admitted inventory")


def build_manifest(config: TrainConfig, *, torch_version: str,
                   cuda: Optional[str]) -> dict:
    """The *pre-fit* run manifest: the declared plan and provenance.

    Written before the optimizer loop, as the charter requires the sampled
    bases and repeats to be declared before fitting. The deterministic draw
    (uniform with replacement, seeded by ``config.seed``) over the train
    windows bounds the workload; ``planned_draws = steps * batch_size`` is the
    declared repeat count and ``max_window`` bounds per-window length. Both
    ``dev_seqids`` and ``eval_every`` are recorded, so two runs that reserve
    different chromosomes or evaluate on different schedules produce distinct
    manifests even at identical result counts (engels-0073).

    Actual attempted/accepted work is attached afterward by ``record_actual``.
    Committed next to the checkpoint (the checkpoint lives on the gagarin
    scratch/return path, not in git)."""
    return {
        "task": "T-human-014",
        "commit": _git_commit(),
        "source": _source_provenance(),
        "param_count": SECTION_35_PARAM_COUNT,
        "hardware": platform.platform(),
        "python": platform.python_version(),
        "torch": torch_version,
        "cuda": cuda,
        "device": config.device,
        # This increment fits the encoder emissions against the fixed-grammar
        # reference chain loss; the learned pooled-decoder scalars are not yet
        # in the loss, so this is an encoder-only profiling fit (engels-0073).
        "scope": "encoder-only-fixed-grammar",
        "hyperparams": {
            "seed": config.seed,
            "steps": config.steps,
            "batch_size": config.batch_size,
            "lr": config.lr,
            "weight_decay": config.weight_decay,
            "grad_clip": config.grad_clip,
            "eval_every": config.eval_every,
            "max_window": config.max_window,
            "loss_kernel": config.loss_kernel,
        },
        "sampling_plan": {
            "draw": "uniform-with-replacement over train windows",
            "seed": config.seed,
            "steps": config.steps,
            "batch_size": config.batch_size,
            "planned_draws": config.steps * config.batch_size,
            "max_window": config.max_window,
        },
        "sources": [
            {"name": s.name, "dev_seqids": list(s.dev_seqids),
             **_source_digests(s)}
            for s in config.sources
        ],
    }


def record_actual(manifest: dict, *, attempted_draws: int, accepted_windows: int,
                  sampled_bases: int, train_windows: int, dev_windows: int,
                  best_dev_nll: Optional[float]) -> dict:
    """Attach the actually-executed work to a pre-fit manifest.

    The charter asks for the planned draw to be declared before fitting and the
    actual attempted/accepted work recorded separately afterward; this keeps
    both in one artifact. ``attempted_draws`` is how many windows the loop drew,
    ``accepted_windows`` how many produced a finite loss (an unusable crop is
    dropped, section 3.6)."""
    out = dict(manifest)
    out["actual"] = {
        "attempted_draws": attempted_draws,
        "accepted_windows": accepted_windows,
        "sampled_bases": sampled_bases,
        "train_windows": train_windows,
        "dev_windows": dev_windows,
        "best_dev_nll": best_dev_nll,
    }
    return out


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


def _window_loss(model, ex, device, dtype, kernel: str = "fast"):
    """Differentiable chain NLL for one window, or ``None`` if the numerator
    support admits no legal path (an unusable crop; section 3.6 says drop it).
    ``kernel`` selects the vectorized delayed-entry forward (``"fast"``, the
    training kernel) or the expanded reference recurrence (``"reference"``)."""
    import torch

    from model.grammar.codes import TABLES

    if kernel == "fast":
        from .fast_loss import chain_nll
    elif kernel == "reference":
        from .torch_loss import chain_nll
    else:
        raise ValueError(f"unknown loss kernel {kernel!r}")

    emissions = _emissions_for(model, ex, device, dtype)
    code = TABLES[ex.table]
    loss = chain_nll(ex.window, emissions, ex.cds_ranges, ex.intron_ranges, code=code)
    if not torch.isfinite(loss):
        return None
    return loss


def train(config: TrainConfig, log=print) -> dict:
    """Fit candidate A's encoder and return the run manifest.

    This is an **encoder-only profiling fit against the fixed-grammar reference
    chain loss**: the loss reads ``model.encoder`` emissions and scores them
    through ``model.grammar``'s default duration/motif model, so the pooled
    decoder scalars (``model.decoder``) have no gradient path and are
    deliberately excluded from the optimizer (engels-0073). Wiring the learned
    pooled duration/motif model into the loss and decoder is the next increment;
    until then this profiles the encoder path and its per-window cost, not the
    full candidate-A objective.

    One gradient step accumulates ``batch_size`` per-window losses (windows vary
    in length, so they are summed rather than tensor-batched; the fast kernel
    will collate). Development NLL is evaluated every ``eval_every`` steps and
    the lowest-NLL checkpoint is kept. The pre-fit manifest (declared plan +
    provenance) is written to ``out_dir/run_manifest.json`` before the loop and
    rewritten with the actual work afterward; the checkpoint is
    ``out_dir/best.pt``.
    """
    import os

    import torch

    from .encoder import CandidateA

    torch.manual_seed(config.seed)
    device = torch.device(config.device)
    dtype = torch.float64  # matches the oracle/parity tests

    examples, stats = _load_all_windows(config)
    # Reject an unusable declared dev split before spending any gradient step.
    validate_dev_reservations(config, stats)
    dev_ids = [i for s in config.sources for i in s.dev_seqids]
    train_ex, dev_ex = split_windows(examples, dev_ids)
    if not train_ex:
        raise ValueError("no training windows after the dev split")
    declared_dev = any(s.dev_seqids for s in config.sources)
    if declared_dev and not dev_ex:  # defensive: validate_dev_reservations covers this
        raise ValueError("declared development reservations retained no windows")
    log(f"loaded {len(examples)} windows: {len(train_ex)} train, {len(dev_ex)} dev")

    model = CandidateA().to(device)
    assert model.num_parameters() == SECTION_35_PARAM_COUNT, model.num_parameters()
    # Encoder-only: the fixed-grammar loss gives model.decoder no gradient.
    opt = torch.optim.Adam(model.encoder.parameters(), lr=config.lr,
                           weight_decay=config.weight_decay)
    gen = torch.Generator().manual_seed(config.seed)

    os.makedirs(config.out_dir, exist_ok=True)
    # Declare the plan and provenance before fitting (charter). Actual work is
    # recorded back into the same file after the loop.
    manifest = build_manifest(
        config, torch_version=torch.__version__,
        cuda=(torch.version.cuda if torch.cuda.is_available() else None))
    manifest_path = os.path.join(config.out_dir, "run_manifest.json")
    with open(manifest_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)

    best_dev = None
    sampled_bases = 0
    attempted_draws = 0
    accepted_windows = 0

    def evaluate() -> Optional[float]:
        if not dev_ex:
            return None
        model.eval()
        total, count = 0.0, 0
        with torch.no_grad():
            for ex in dev_ex:
                loss = _window_loss(model, ex, device, dtype, config.loss_kernel)
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
            attempted_draws += 1
            loss = _window_loss(model, ex, device, dtype, config.loss_kernel)
            if loss is None:
                continue
            (loss / config.batch_size).backward()
            batch_loss += float(loss.detach())
            sampled_bases += ex.n
            accepted_windows += 1
            used += 1
        if used:
            torch.nn.utils.clip_grad_norm_(model.encoder.parameters(), config.grad_clip)
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

    manifest = record_actual(
        manifest, attempted_draws=attempted_draws,
        accepted_windows=accepted_windows, sampled_bases=sampled_bases,
        train_windows=len(train_ex), dev_windows=len(dev_ex), best_dev_nll=best_dev)
    with open(manifest_path, "w", encoding="utf-8") as fh:
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
            checkpoint: Optional[str], json_out: Optional[str] = None,
            decoder: str = "tensor", decode_batch: int = 16, log=print) -> dict:
    """Time preprocessing, encoder and decode/traceback over one sequence.

    Returns a dict of measured stage times, oriented bases, and per-Mb rates,
    plus peak host RSS and device memory. Each stage records **both** process
    CPU seconds and synchronized elapsed wall seconds: on ``cuda`` the device
    runs asynchronously, so ``time.process_time`` does not see kernel execution;
    the GPU budget (``gpu_s_per_mb``) is taken from the CUDA-synchronized wall
    clock, the CPU budget (``cpu_s_per_mb``) from process CPU (engels-0073).

    This is an **annotation-selected window profile**, not an end-to-end
    chromosome row: ``_load_all_windows`` runs before timing, only clean
    admitted gene windows enter the denominator, and predictions are discarded
    (``outputs_discarded``). It profiles the encoder/decoder cost per admitted
    base; it does not include full-sequence, both-strand, I/O, output-writing or
    multi-worker accounting, so it cannot fill a chromosome sensitivity/budget
    row on its own. The decoder is the fixed-grammar delayed-entry Viterbi,
    consistent with the encoder-only training scope; the learned pooled
    decoder is not yet wired. ``decoder`` selects the tensor max-product scan
    (``"tensor"``, :mod:`model.a.fast_viterbi`, what ships) or the
    pure-Python semantics check (``"python"``, the 66 CPU-s/Mb stage of
    a-pilot.md section 3.2, kept for parity runs). The tensor scan decodes
    the windows in length-sorted batches of ``decode_batch`` (per-window
    encoder timing is unchanged; the decode stage is timed as a whole), which
    is how a chromosome's windows would be decoded in either regime.
    """
    import torch

    from .encoder import CandidateA
    from .fast_viterbi import viterbi_windows
    from model.grammar.codes import TABLES
    from model.grammar.delayed import DelayedEntryDecoder

    if decoder not in ("tensor", "python"):
        raise ValueError("decoder must be 'tensor' or 'python'")

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

    def _sync():
        if device.type == "cuda":
            torch.cuda.synchronize(device)

    pre_cpu = enc_cpu = dec_cpu = 0.0
    pre_wall = enc_wall = dec_wall = 0.0
    bases = 0
    pending_windows, pending_emissions, pending_codes = [], [], []
    with torch.no_grad():
        for ex in examples:
            padded, available = pad_window(ex.window)
            c0, w0 = time.process_time(), time.perf_counter()
            feats = encode_sequence(padded, available=available).unsqueeze(0).to(device)
            _sync()
            pre_cpu += time.process_time() - c0
            pre_wall += time.perf_counter() - w0

            c0, w0 = time.process_time(), time.perf_counter()
            emissions = model.encoder(feats)[0][:, : ex.n].to(dtype)
            _sync()
            enc_cpu += time.process_time() - c0
            enc_wall += time.perf_counter() - w0

            if decoder == "tensor":
                pending_windows.append(ex.window)
                pending_emissions.append(emissions)
                pending_codes.append(TABLES[ex.table])
            else:
                c0, w0 = time.process_time(), time.perf_counter()
                sc = _scores_from_emissions(emissions)
                DelayedEntryDecoder(TABLES[ex.table]).viterbi(ex.window, sc)
                dec_cpu += time.process_time() - c0
                dec_wall += time.perf_counter() - w0
            bases += ex.n
        if decoder == "tensor":
            c0, w0 = time.process_time(), time.perf_counter()
            viterbi_windows(pending_windows, pending_emissions, codes=pending_codes,
                            batch_size=decode_batch)
            _sync()
            dec_cpu += time.process_time() - c0
            dec_wall += time.perf_counter() - w0

    mb = bases / 1e6
    device_mem_gib = (torch.cuda.max_memory_allocated(device) / 2**30
                     if device.type == "cuda" else None)
    cpu_total = pre_cpu + enc_cpu + dec_cpu
    wall_total = pre_wall + enc_wall + dec_wall
    # The GPU budget is elapsed device time (CUDA-synchronized wall), not CPU.
    gpu_s_per_mb = wall_total / mb if (device.type == "cuda" and mb) else None
    row = {
        "species": species, "seqid": seqid, "windows": len(examples),
        "oriented_bases": bases, "device": config.device,
        "profile": "annotation-selected-windows", "outputs_discarded": True,
        "decoder": decoder, "decode_batch": decode_batch if decoder == "tensor" else 1,
        "commit": _git_commit(), "source": _source_provenance(),
        "preprocess_cpu_s": pre_cpu, "encoder_cpu_s": enc_cpu, "decode_cpu_s": dec_cpu,
        "preprocess_wall_s": pre_wall, "encoder_wall_s": enc_wall,
        "decode_wall_s": dec_wall,
        "cpu_s_per_mb": cpu_total / mb if mb else None,
        "wall_s_per_mb": wall_total / mb if mb else None,
        "gpu_s_per_mb": gpu_s_per_mb,
        # ``ru_maxrss`` is KiB on Linux; both peaks are reported in GiB (2^30
        # bytes) and named so (engels-0080 P3), the convention a-pilot.md uses.
        "peak_host_rss_gib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20,
        "peak_device_mem_gib": device_mem_gib,
    }
    log(json.dumps(row, indent=2, sort_keys=True))
    if json_out:
        # A complete machine-readable artifact, independent of stdout capture
        # length or interleaving with the /usr/bin/time report.
        with open(json_out, "w") as fh:
            json.dump(row, fh, indent=2, sort_keys=True)
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
    pm.add_argument("--json-out", default=None,
                    help="also write the complete measurement row to this file")
    pm.add_argument("--decoder", default="tensor", choices=("tensor", "python"),
                    help="Viterbi implementation to time (default: tensor scan)")
    pm.add_argument("--decode-batch", type=int, default=16,
                    help="windows per length-sorted decoding batch (tensor decoder)")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    config = TrainConfig.from_json(args.config)
    if args.cmd == "train":
        train(config)
    elif args.cmd == "measure":
        measure(config, args.species, args.seqid, args.checkpoint,
                json_out=args.json_out, decoder=args.decoder,
                decode_batch=args.decode_batch)
    return 0


if __name__ == "__main__":
    sys.exit(main())
