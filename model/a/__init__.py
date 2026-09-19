"""Candidate A: compact DNA encoder and duration-mixture CRF (proposal section 3).

This package holds the *neural* half of candidate A: the eight-channel
featurizer (`features`) and the encoder that maps a DNA chunk to the eleven
per-base emission channels the `model.grammar` decoders consume, plus the
54 pooled-decoder scalars. The decoders themselves (reference and delayed
entry) already live in `model.grammar` from T-human-013.

The parameter inventory is fixed by proposal section 3.5 at exactly
455,841 trainable scalars. `SECTION_35_PARAM_COUNT` reproduces that number
from the formulas so a test can guard both the spec and the built module.

Fitting and end-to-end measurement (T-human-014's remaining deliverables)
run on gagarin under the Phase 4 budget; nothing here trains a model.
"""
from .inventory import (
    DEPENDENCY_RADIUS,
    SECTION_35_PARAM_COUNT,
    section_35_param_breakdown,
)
from .loss import chain_nll, numerator_scores

# The data loader (`model.a.dataset`) pulls in `model.labels`, whose audit loads
# `benchmark/score.py` by path: it is a checkout-time tool, not part of an
# installed wheel's import graph. Import it lazily (PEP 562) so `import model.a`
# -- and the torch-free parameter-count guard, the encoder, and the loss -- do
# not require `model.labels` or `benchmark/` to be present.
_LAZY = {"LoaderStats", "SourceMismatch", "WindowExample", "iter_windows", "verify_source"}


def __getattr__(name):
    if name in _LAZY:
        from . import dataset
        return getattr(dataset, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "SECTION_35_PARAM_COUNT",
    "section_35_param_breakdown",
    "DEPENDENCY_RADIUS",
    "chain_nll",
    "numerator_scores",
    "iter_windows",
    "verify_source",
    "WindowExample",
    "LoaderStats",
    "SourceMismatch",
]

# `CandidateA`, `DNAEncoder`, `DecoderParams` are torch modules; import them
# lazily so this package (and the stdlib arithmetic guard) works without torch.
try:  # pragma: no cover - exercised only where torch is installed
    from .encoder import CandidateA, DNAEncoder, DecoderParams  # noqa: F401

    __all__ += ["CandidateA", "DNAEncoder", "DecoderParams"]
except ModuleNotFoundError:  # pragma: no cover
    pass
