"""Genetic-code tables the grammar consumes (proposal section 3.1).

Only the stop set and the allowed initiators matter to the decoder: S
states accept a codon as an initiator, E states terminate the chain on a
stop. Tables follow the NCBI genetic codes page
(https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). Alternative
initiators are listed per table but not enabled unless a declaration asks
for them, since a permissive initiator set only widens the grammar.
"""
from dataclasses import dataclass, field

# IUPAC ambiguity codes -> permitted bases (proposal 3.1: a normalized prior
# over the permitted bases, uniform to begin with).
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT", "X": "ACGT",
}


def permitted_bases(x):
    """Bases an observed symbol may stand for; unknown symbols behave like N."""
    return IUPAC.get(x.upper(), "ACGT")


@dataclass(frozen=True)
class GeneticCode:
    table: int
    name: str
    stops: frozenset
    initiators: frozenset            # enabled initiators
    alternative_initiators: frozenset = field(default_factory=frozenset)

    def with_alternative_initiators(self):
        return GeneticCode(self.table, self.name, self.stops,
                           self.initiators | self.alternative_initiators,
                           self.alternative_initiators)

    def is_stop(self, codon):
        return codon in self.stops

    def is_initiator(self, codon):
        return codon in self.initiators

    def initiator_prefixes(self):
        """All proper prefixes (length 1 and 2) of enabled initiators."""
        return frozenset(c[:k] for c in self.initiators for k in (1, 2))


TABLES = {
    1: GeneticCode(1, "Standard", frozenset({"TAA", "TAG", "TGA"}), frozenset({"ATG"}),
                   frozenset({"TTG", "CTG"})),
    # Ciliate, dasycladacean and hexamita nuclear: TAA and TAG code Gln; only
    # TGA stops (proposal 3.1, acceptance case 4).
    6: GeneticCode(6, "Ciliate nuclear", frozenset({"TGA"}), frozenset({"ATG"})),
}
