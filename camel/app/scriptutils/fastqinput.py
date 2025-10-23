from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.snakemake import snakemakeutils


@dataclass(frozen=True)
class FastqInput:
    """
    Dataclass used to store FASTQ input/output in a standardized format.
    """
    read_type: str
    pe: Optional[list[ToolIOFile]] = None
    se: Optional[list[ToolIOFile]] = None
    se_fwd: Optional[list[ToolIOFile]] = None
    se_rev: Optional[list[ToolIOFile]] = None
    is_trimmed: bool = False
    is_pe: bool = True

    @property
    def is_paired(self) -> bool:
        """
        Returns true if the read input is paired.
        :return: True if paired, False otherwise
        """
        return self.pe is not None

    def to_fq_dict(self) -> dict[str, list[ToolIOFile]]:
        """
        Converts the FASTQ input to an input dictionary for the workflow.
        :return: FASTQ dictionary
        """
        if self.read_type == 'hybrid':
            fq_dict = {'PE': self.pe, 'SE': self.se}
            if self.se_fwd is not None:
                fq_dict['SE_FWD'] = self.se_fwd
            if self.se_rev is not None:
                fq_dict['SE_REV'] = self.se_rev
            return fq_dict
        elif self.is_pe:
            fq_dict = {'PE': self.pe}
            if self.se_fwd is not None:
                fq_dict['SE_FWD'] = self.se_fwd
            if self.se_rev is not None:
                fq_dict['SE_REV'] = self.se_rev
            return fq_dict
        else:
            return {'SE': self.se}

    @staticmethod
    def from_fq_dict(io: Path, read_type: str) -> 'FastqInput':
        """
        Creates a FastqInput from an IO object.
        :param io: IO object
        :param read_type: Read type
        :return: FastqInput
        """
        fq_dict = snakemakeutils.load_object(io)
        return FastqInput(
            read_type,
            pe=fq_dict.get('PE'),
            se_fwd=fq_dict.get('SE_FWD'),
            se_rev=fq_dict.get('SE_REV'),
            se=fq_dict.get('SE'),
            is_pe=fq_dict.get('SE') is None
        )
