import re
from Bio.PDB import Polypeptide
# noinspection PyProtectedMember
from vcf.model import _Record as VcfRecord

from camel.app.toolkits.csq.csqutils import BCSQInfo
from camel.app.toolkits.csq.mutations.basemutation import BaseMutation


class AminoAcidMutation(BaseMutation):
    """
    This class represents a mutation that leads to a missense or a synonymous mutation.
    """

    def __init__(self, record: VcfRecord, locus: str, pos_rel: int, aa1: str, aa2: str) -> None:
        """
        Initializes an amino acid mutation.
        :param record: VCF record
        :param locus: Locus
        :param pos_rel: Relative position
        :param aa1: Reference amino acid
        :param aa2: Alternate amino acid
        """
        super().__init__(record)
        self.locus = locus
        self.pos_rel = pos_rel
        self.aa1 = aa1
        self.aa2 = aa2

    @staticmethod
    def parse(record: VcfRecord, info: BCSQInfo) -> 'AminoAcidMutation':
        """
        Parses an amino acid mutation.
        :param record: VCF record
        :param info: BCSQ information
        :return: Parsed mutation
        """
        parts = info.raw_str.split('|')
        # Check for AA substitution
        m = re.match('(\\d+)(\\w)>\\d+(\\w)', parts[5])
        if m:
            return AminoAcidMutation(record, parts[2], int(m.group(1)), m.group(2), m.group(3))

        # Check for synonymous
        m = re.match('^(\\d+)(\\w)$', parts[5])
        if m:
            return AminoAcidMutation(record, parts[2], int(m.group(1)), m.group(2), m.group(2))

        raise ValueError(f"Cannot parse {AminoAcidMutation.__class__.__name__} mutation from '{info}'")

    @property
    def short_notation(self) -> str:
        """
        Returns the short notation of this mutation.
        :return: Short notation
        """
        return ''.join([self.aa1_code_1, str(self.pos_rel), self.aa2_code_1])

    @property
    def long_notation(self) -> str:
        """
        Returns the long notation of this mutation.
        :return: Long notation
        """
        return '-'.join([self.aa1_code_3, str(self.pos_rel), self.aa2_code_3])

    @property
    def aa1_code_1(self) -> str:
        """
        Returns the first amino acid (1 letter notation)
        :return: First amino acid
        """
        return self.aa1

    @property
    def aa1_code_3(self) -> str:
        """
        Returns the first amino acid (3 letter notation)
        :return: First amino acid
        """
        return Polypeptide.one_to_three(self.aa1).title()

    @property
    def aa2_code_1(self) -> str:
        """
        Returns the second amino acid (1 letter notation)
        :return: Second amino acid
        """
        return self.aa2

    @property
    def aa2_code_3(self) -> str:
        """
        Returns the second amino acid (3 letter notation)
        :return: Second amino acid
        """
        return Polypeptide.one_to_three(self.aa2).title()

    def is_synonymous(self) -> bool:
        """
        Returns true if the mutation is synonymous (same amino-acid)
        :return: True if synonymous, False otherwise
        """
        return self.aa1 == self.aa2

    @property
    def key(self):
        """
        Returns the key of the mutation, which can be used to compare against the database mutations.
        :return: Key
        """
        return 'AA', self.aa1, self.pos_rel, self.aa2
