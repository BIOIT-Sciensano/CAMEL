import re
# noinspection PyProtectedMember
from vcf.model import _Record as VcfRecord

from camel.app.toolkits.csq.csqutils import BCSQInfo
from camel.app.toolkits.csq.mutations.basemutation import BaseMutation


class StopMutation(BaseMutation):
    """
    This class contains a mutation that leads to a stop codon inside the coding region.
    """

    def __init__(self, record: VcfRecord, locus: str, pos_rel: int) -> None:
        """
        Initializes a stop codon mutation.
        :param record: VCF record
        :param locus: Locus
        :param pos_rel: Relative position
        """
        super().__init__(record)
        self.locus = locus
        self.pos_rel = pos_rel

    @staticmethod
    def parse(record: VcfRecord, info: BCSQInfo) -> 'StopMutation':
        """
        Parses a stop mutation.
        :param record: VCF record
        :param info: BCSQ information
        :return: Parsed mutation
        """
        parts = info.raw_str.split('|')
        m = re.match('^(\\d+)\\w', parts[5])
        if m:
            return StopMutation(record, parts[2], int(m.group(1)))
        raise ValueError(f"Cannot parse {StopMutation.__class__.__name__} mutation from '{info}'")

    @property
    def short_notation(self) -> str:
        """
        Returns the short notation of this mutation.
        :return: Short notation
        """
        return 'stop'

    @property
    def long_notation(self) -> str:
        """
        Returns the long notation of this mutation.
        :return: Long notation
        """
        return f'Stop {self.pos_rel}'

    @property
    def key(self):
        """
        Returns the key of the mutation, which can be used to compare against the database mutations.
        :return: Key
        """
        return 'STOP',
