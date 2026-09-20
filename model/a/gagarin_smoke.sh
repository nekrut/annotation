#!/usr/bin/env bash
# Bounded verification + profiling run for candidate A (T-human-014).
#
# Purpose (NOT the full pilot fit): on a real torch host,
#   1. run the full candidate-A torch test suite -- the input-validation,
#      parity and autograd tests that skip on the torch-free dev host -- to
#      verify the tensor path,
#   2. fetch and checksum-verify the smallest train species (yeast) and run a
#      short `train` + `measure` to profile the reference-recurrence per-window
#      cost of the differentiable chain loss.
#
# gagarin runs exactly this script. It assumes a Python 3.11 interpreter with
# torch importable (the cluster's Phase 4 env) and network access to NCBI's
# genomes/all for the pinned yeast source. Everything it writes goes under
# $SCRATCH; nothing is committed here.
#
# Usage: bash model/a/gagarin_smoke.sh [SCRATCH_DIR]
set -euo pipefail

REPO="$(git rev-parse --show-toplevel)"
cd "$REPO"
SCRATCH="${1:-${SCRATCH:-$(mktemp -d)}}"
SOURCES="$SCRATCH/sources"
OUT="$SCRATCH/run"
mkdir -p "$SOURCES" "$OUT"
export PYTHONPATH="$REPO:${PYTHONPATH:-}"

echo "== commit =="
git rev-parse HEAD

echo "== torch env =="
python3 -c "import torch, platform; print('torch', torch.__version__, 'cuda', torch.version.cuda, torch.cuda.is_available()); print(platform.platform())"

echo "== 1. full candidate-A torch test suite =="
python3 -m unittest -v \
  tests.test_a_encoder tests.test_a_loss tests.test_a_torch_loss \
  tests.test_a_dataset tests.test_a_train

echo "== 2a. fetch + checksum-verify the smallest train species (yeast) =="
python3 -m model.a.coverage --fetch --species Saccharomyces_cerevisiae \
  --sources "$SOURCES" --out "$OUT/coverage_yeast.tsv"
cat "$OUT/coverage_yeast.tsv"

echo "== 2b. write the smoke TrainConfig =="
DEVICE="cpu"; python3 -c "import torch,sys; sys.exit(0 if torch.cuda.is_available() else 1)" && DEVICE="cuda"
python3 - "$SOURCES" "$OUT" "$DEVICE" <<'PY'
import json, os, sys
from model.a import coverage as C
from model.a.coverage import species_sources
sources_dir, out, device = sys.argv[1], sys.argv[2], sys.argv[3]
here = os.path.dirname(os.path.abspath(C.__file__))
manifests_dir = os.path.normpath(os.path.join(here, "..", "labels", "manifests"))
srcs = [s for s in species_sources(manifests_dir, sources_dir) if s[0] == "Saccharomyces_cerevisiae"]
assert srcs, f"no yeast summary under {manifests_dir}"
sp, summary, gff, fasta = srcs[0]
cfg = {
    "sources": [{"name": sp, "summary": summary, "gff": gff, "fasta": fasta, "dev_seqids": []}],
    "out_dir": out, "seed": 0, "steps": 20, "batch_size": 4, "lr": 3e-4,
    "eval_every": 10, "max_window": 12288, "device": device,
}
with open(os.path.join(out, "smoke_config.json"), "w") as fh:
    json.dump(cfg, fh, indent=2)
print("wrote", os.path.join(out, "smoke_config.json"), "device", device)
PY

echo "== smoke config (returned so the split/plan is reconstructable) =="
cat "$OUT/smoke_config.json"

echo "== 2c. short train (20 steps) — profiles the reference-recurrence loss cost =="
/usr/bin/time -v python3 -m model.a.train train --config "$OUT/smoke_config.json" 2>&1 | tail -40
echo "== run manifest (declared plan + actual work; carries the split/schedule) =="
cat "$OUT/run_manifest.json"

echo "== 2d. measure preprocessing / encoder / decode on yeast =="
/usr/bin/time -v python3 -m model.a.train measure --config "$OUT/smoke_config.json" \
  --species Saccharomyces_cerevisiae --checkpoint "$OUT/best.pt" 2>&1 | tail -40

echo "== done; artifacts under $OUT =="
