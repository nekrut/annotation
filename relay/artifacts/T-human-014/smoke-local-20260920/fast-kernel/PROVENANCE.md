# Provenance of the fast-kernel records in this directory

Written 2026-09-20 in response to engels-0080 (P2): the raw manifests below are
preserved unchanged; this file reconciles what they record with what ran.

| record | manifest `commit` | tree state when run | source that ran | exact pin? |
|---|---|---|---|---|
| `run_manifest.json`, `train.out`, `train_time.txt` (fast fit, default threads) | `751700b` | dirty: `model/a/fast_loss.py` (untracked), `model/a/train.py` (`loss_kernel` field + fast path), `tests/test_a_fast_loss.py` | the working tree later committed as `28c5339` in the same session; no edit between run and commit is recorded, but no hash of the executed tree was taken | **no** |
| `one-thread/` (fast fit, one core) | `751700b` | same | same | **no** |
| `profile_window_fast.out` | n/a (script prints no SHA) | same | same | **no** |
| `pinned-129dcbe/` (fast fit, one core, re-run) | `129dcbe` | clean (`source_dirty: false`) | `129dcbe` exactly; `source_sha256` recorded in the manifest | **yes** |

`751700b` has neither `model/a/fast_loss.py` nor the `loss_kernel` config field,
so checking it out cannot reproduce the first three records; they are kept as
the observed measurements with this qualification. The re-run in
`pinned-129dcbe/` reproduces the one-core number at a committed source
(23.27 s wall / 23.27 CPU-s, 914,560 KiB = 0.87 GiB, NLL 0.0001 / 0.0006 at
steps 10 / 20, identical to `one-thread/`). `129dcbe` differs from `28c5339`
only by the cache-key fix, the `dtype=torch.long` symbol tensors, the manifest
provenance fields and the GiB key names; none touches the ATG-only path the
smoke fit exercises.

From `129dcbe` on, `model/a/train.py` writes `source.source_sha256` (over
`model/a/*.py` + `model/grammar/*.py`), `source.source_dirty` and
`source.source_dirty_files` into every manifest next to `commit`.

Paths in `pinned-129dcbe/config.json` are redacted to `<checkout>` /
`<scratch>` placeholders; the data are the MD5-pinned yeast sources of
`run_manifest.json`.
