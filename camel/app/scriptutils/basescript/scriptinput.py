import dataclasses
from pathlib import Path
from typing import Any

from camel.app.core.utils import fastautils
from camel.app.scriptutils import model


@dataclasses.dataclass(frozen=True)
class ScriptInput(model.BaseInput):
    """
    Contains the input for a script or pipeline.
    """

    type_: model.InputType
    sample_name: str | None = None
    fastq_se: Path | None = None
    fastq_se_name: str | None = None
    fastq_pe: tuple[Path, Path] | None = None
    fastq_pe_names: tuple[str, str] | None = None
    fasta: Path | None = None
    fasta_name: str | None = None
    vcf_unfiltered: Path | None = None

    @property
    def input_str(self) -> str:
        """
        Returns the input file string for this input.
        :return: Input file string
        """
        if self.type_ == model.InputType.FASTA:
            return self.fasta_name if self.fasta_name else self.fasta.name
        elif self.type_ == model.InputType.FASTA_WITH_VCF:
            name_fasta = self.fasta_name if self.fasta_name else self.fasta.name
            name_vcf = self.vcf_unfiltered.name
            return f'{name_fasta}, {name_vcf}'
        elif self.type_ == model.InputType.ONT:
            return self.fastq_se_name if self.fastq_se_name else self.fastq_se.name
        elif self.type_ == model.InputType.ILLUMINA:
            parts = self.fastq_pe_names if self.fastq_pe_names else [p.name for p in self.fastq_pe]
            return ', '.join(parts)
        elif self.type_ == model.InputType.HYBRID:
            parts = list(self.fastq_pe_names) if self.fastq_pe_names else [p.name for p in self.fastq_pe]
            parts.append(self.fastq_se_name if self.fastq_se_name else self.fastq_se.name)
            return ', '.join(parts)
        raise ValueError(f"Unknown input type: {self.type_}")

    @property
    def name(self) -> str:
        """
        Returns the dataset name.
        :return: Dataset name
        """
        if self.sample_name is not None:
            return self.sample_name
        return "NA"

    def to_dict(self) -> dict[str, Any]:
        """
        Converts the script input to a dictionary.
        :return: Dictionary with the input files
        """
        data: dict[str, Any] = {"sample_name": self.sample_name, "type": self.type_.value}

        # FASTA
        if self.fasta is not None:
            data["fasta"] = ([{
                "name": self.fasta_name if self.fasta_name is not None else self.fasta.name,
                "path": str(self.fasta)
            }])

        # FASTQ SE
        if self.fastq_se is not None:
            data["fastq_se"] = [{
                "name": self.fastq_se_name if self.fastq_se_name is None else self.fastq_se_name,
                "path": str(self.fastq_se)}
            ]

        # FASTQ PE
        if self.fastq_pe is not None:
            data["fastq_pe"] = [{"name": p.name, "path": str(p)} for p in self.fastq_pe]

        # VCF
        if self.vcf_unfiltered is not None:
            data["vcf"] = [{"name": self.vcf_unfiltered.name, "path": str(self.vcf_unfiltered)}]

        return data

    @staticmethod
    def from_dict(data_in: dict[str, Any]) -> 'ScriptInput':
        """
        Restores the script input from a dictionary.
        :param data_in: Dictionary with the input files
        :return: Script input
        """
        if 'fasta' in data_in:
            if not len(data_in['fasta']) == 1:
                raise ValueError('Only one FASTA file is supported.')
        if 'fastq_se' in data_in:
            if not len(data_in['fastq_se']) == 1:
                raise ValueError('Only one single-end FASTQ file is supported.')
        if 'fastq_pe' in data_in:
            if not len(data_in['fastq_pe']) == 2:
                raise ValueError('Only two paired-end FASTQ files are supported.')
        if 'vcf' in data_in:
            if not len(data_in['vcf']) == 1:
                raise ValueError('Only one VCF file is supported.')

        # Returns the output
        return ScriptInput(
            type_=model.InputType(data_in['type']),
            sample_name=data_in['sample_name'],
            # FASTA
            fasta=Path(data_in['fasta'][0]['path']) if 'fasta' in data_in else None,
            fasta_name=data_in['fasta'][0]['name'] if 'fasta' in data_in else None,

            # FASTQ SE
            fastq_se=Path(data_in['fastq_se'][0]['path']) if 'fastq_se' in data_in else None,
            fastq_se_name=data_in['fastq_se'][0]['name'] if 'fastq_se' in data_in else None,

            # FASTQ PE
            fastq_pe=(
                Path(data_in['fastq_pe'][0]['path']),
                Path(data_in['fastq_pe'][1]['path'])
            ) if 'fastq_pe' in data_in else None,
            fastq_pe_names=(
                data_in['fastq_pe'][0]['name'],
                data_in['fastq_pe'][1]['name']
            ) if 'fastq_pe' in data_in else None,

            # VCF file
            vcf_unfiltered=Path(data_in['vcf'][0]['path']) if 'vcf' in data_in else None,
        )

    def validate(self) -> None:
        """
        Validates if the script input is valid.
        """
        if self.type_ == model.InputType.FASTA:
            if fastautils.has_duplicates(self.fasta):
                raise ValueError('Input FASTA file has duplicate sequence IDs.')
