from typing import Optional

# noinspection PyProtectedMember
from vcf.model import _Record as VcfRecord

from camel.app.components.csq.mutations.basemutation import BaseMutation


class NucleotideMutation(BaseMutation):
    """
    This class contains a mutation in a non coding region.
    """

    def __init__(self, record: VcfRecord, pos_rel: Optional[int], locus_tabix: str) -> None:
        """
        Initializes a nucleotide mutation.
        :param record: VCF record
        :param pos_rel: Relative position
        :param locus_tabix: Tabix locus name
        """
        super().__init__(record)
        self.locus_tabix = locus_tabix
        self.pos_rel = pos_rel

    @staticmethod
    def parse(record: VcfRecord, tabix_by_pos: dict[int, dict]) -> 'NucleotideMutation':
        """
        Parses a nucleotide mutation.
        :param record: VCF record
        :param tabix_by_pos: TABIX annotation by position
        :return: Parsed mutation
        """
        tabix_annot = tabix_by_pos[record.POS]
        pos_rel = int(tabix_annot['mutation_pos']) if tabix_annot['mutation_pos'] != '.' else None
        return NucleotideMutation(record, pos_rel, tabix_annot['locus'])

    @property
    def short_notation(self) -> str:
        """
        Returns the short notation of this mutation.
        :return: Short notation
        """
        return ''.join([
            self.ref,
            '({})'.format(str(self.pos_rel) if self.pos_rel is not None else '.'),
            ', '.join([n for n in self.alt])
        ])

    @property
    def long_notation(self) -> str:
        """
        Returns the long notation of this mutation.
        :return: Long notation
        """
        return ''.join([
            self.ref.lower(),
            '({})'.format(str(self.pos_rel) if self.pos_rel is not None else '.'),
            ', '.join([n.lower() for n in self.alt])
        ])

    @property
    def key(self) -> tuple[str, str, int, str]:
        """
        Returns the mutation key.
        """
        return 'NUCL', self.ref.lower(), self.pos_rel, str(self.alt[0]).lower()
