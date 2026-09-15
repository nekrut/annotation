"""Candidate A's coding grammar: state inventory, genetic codes and decoders.

`reference.py` is the expanded, slow recurrence that proposal section 3.2
names as the specification oracle; the delayed-entry decoder is checked
against it. Standard library only.
"""
from .codes import GeneticCode, TABLES, permitted_bases
from .scores import Scores
from .reference import ReferenceDecoder, DurationMixture, Chain, Segment

__all__ = ["GeneticCode", "TABLES", "permitted_bases", "Scores",
           "ReferenceDecoder", "DurationMixture", "Chain", "Segment"]
