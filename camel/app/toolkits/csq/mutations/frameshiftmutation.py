import re
# noinspection PyProtectedMember
from vcf.model import _Record as VcfRecord

from camel.app.toolkits.csq.csqutils import BCSQInfo
from camel.app.toolkits.csq.mutations.basemutation import BaseMutation


class FrameshiftMutation(BaseMutation):
    """
    This class represents a mutation that leads to a frameshift.
    """

    def __init__(self, record: VcfRecord, locus: str, pos_rel: int) -> None:
        """
        Initializes a frameshift mutation.
        :param record: VCF record
        :param locus: Locus
        :param pos_rel: Relative position
        """
        super().__init__(record)
        self.locus = locus
        self.pos_rel = pos_rel

    @staticmethod
    def parse(record: VcfRecord, info: BCSQInfo) -> 'FrameshiftMutation':
        """
        Parses a frameshift mutation.
        :param record: VCF record
        :param info: BCSQ information
        :return: Parsed mutation
        """
        parts = info.raw_str.split('|')
        m = re.match('^(\\d+)\\w', parts[5])
        if m:
            return FrameshiftMutation(record, parts[2], int(m.group(1)))
        raise ValueError(f"Cannot parse {FrameshiftMutation.__class__.__name__} mutation from '{info}'")

    @property
    def short_notation(self) -> str:
        """
        Returns the short notation of this mutation.
        :return: Short notation
        """
        return 'frameshift'

    @property
    def long_notation(self) -> str:
        """
        Returns the long notation of this mutation.
        :return: Long notation
        """
        return f'Frameshift {self.pos_rel}'

    @property
    def key(self):
        """
        Returns the key of the mutation, which can be used to compare against the database mutations.
        :return: Key
        """
        return 'FRAMESHIFT',
