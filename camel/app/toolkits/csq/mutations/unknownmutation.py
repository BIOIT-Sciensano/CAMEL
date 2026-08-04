# noinspection PyProtectedMember
from vcf.model import _Record as VcfRecord

from camel.app.toolkits.csq.csqutils import BCSQInfo
from camel.app.toolkits.csq.mutations.basemutation import BaseMutation


class UnknownMutation(BaseMutation):
    """
    This class contains a mutation not supported by the parser.
    """

    def __init__(self, record: VcfRecord, type_: str) -> None:
        """
        Initializes this mutation.
        :param record: VCF record
        :param type_: Mutation type
        """
        super().__init__(record)
        self.type_ = type_

    @staticmethod
    def parse(record: VcfRecord, info: BCSQInfo) -> 'UnknownMutation':
        """
        Parses an amino acid mutation.
        :param record: VCF record
        :param info: BCSQ information
        :return: Parsed mutation
        """
        return UnknownMutation(record, info.type_)

    @property
    def short_notation(self) -> str:
        """
        Returns the short notation of this mutation.
        :return: Short notation
        """
        return ''.join([self.type_, str(self.start)])

    @property
    def long_notation(self) -> str:
        """
        Returns the long notation of this mutation.
        :return: Long notation
        """
        return f'{self.type_}@{self.start}'
