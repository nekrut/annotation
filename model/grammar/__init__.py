"""Candidate A's coding grammar: state inventory, genetic codes and decoders.

`reference.py` is the expanded, slow recurrence that proposal section 3.2
names as the specification oracle; `delayed.py` is the delayed-entry
recurrence checked against it. Standard library only.
"""
from .codes import GeneticCode, TABLES, permitted_bases
from .scores import Scores
from .reference import ReferenceDecoder, DurationMixture, EdgePrior, Chain, Segment
from .delayed import DelayedEntryDecoder
from .strand import reverse_complement, genomic_features, codon_features, gff3_rows, Feature

__all__ = ["GeneticCode", "TABLES", "permitted_bases", "Scores",
           "ReferenceDecoder", "DelayedEntryDecoder", "DurationMixture",
           "EdgePrior", "Chain", "Segment",
           "reverse_complement", "genomic_features", "codon_features", "gff3_rows", "Feature"]
