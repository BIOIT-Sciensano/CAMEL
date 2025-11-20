import abc
import dataclasses
from enum import Enum


class InputType(Enum):
    """
    Enumerator for input types.
    """

    FASTA = "fasta"
    FASTA_WITH_VCF = "fasta_with_vcf"
    HYBRID = "hybrid"
    ILLUMINA = "illumina"
    ONT = "ont"


@dataclasses.dataclass(frozen=True)
class BaseInput(metaclass=abc.ABCMeta):
    """
    Base class for tool input.
    """

    def validate(self) -> bool:
        """
        Checks if this script input is valid.
        """
        return True

    def name(self) -> str:
        """
        Returns the dataset name.
        """
        return 'NA'


@dataclasses.dataclass(frozen=True)
class BaseOutput(metaclass=abc.ABCMeta):
    """
    Base class for tool output.
    """

    pass


@dataclasses.dataclass(frozen=True)
class BaseOptions:
    """
    Base class for tool options.
    """

    pass
