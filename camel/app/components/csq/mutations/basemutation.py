import abc
from abc import ABCMeta
# noinspection PyProtectedMember
from vcf.model import _Record as VcfRecord


class BaseMutation(metaclass=ABCMeta):
    """
    This class is the base class for mutations supported by csq parser.
    """

    def __init__(self, record: VcfRecord) -> None:
        """
        Initializes the base mutation class.
        """
        self.chrom = record.CHROM
        self.start = record.start
        self.end = record.end
        self.ref = str(record.REF)
        self.alt = [str(x) for x in record.ALT]

    @property
    @abc.abstractmethod
    def long_notation(self) -> str:
        """
        Returns the long notation of this mutation.
        :return: Long notation
        """
        pass

    @property
    @abc.abstractmethod
    def short_notation(self) -> str:
        """
        Returns the short notation of this mutation.
        :return: Short notation
        """
        pass

    @property
    def key(self) -> tuple:
        """
        Returns the key of the mutation, which can be used to compare against the database mutations.
        :return: Key
        """
        return 'key',

    def matches_record(self, record: VcfRecord) -> bool:
        """
        Return true if the mutation matches the given VCF record.
        :return: True if match, False otherwise
        """
        if self.start != record.start:
            return False
        if self.end != record.end:
            return False
        if self.chrom != record.CHROM:
            return False
        if self.alt != [str(x) for x in record.ALT]:
            return False
        return True
