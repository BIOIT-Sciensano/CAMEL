import dataclasses
import logging
import tempfile
from pathlib import Path

from Bio import SeqIO

from camel.app.camel import Camel
from camel.app.command.command import Command


@dataclasses.dataclass
class ApplyVariantsOutput:
    """
    Holder for the output of the 'apply variants' script.
    """
    path_fasta: Path
    informs: list[dict]


class ApplyVariants:
    """
    Wrapper that applies the variants from a VCF file to a FASTA file.
    """

    def __init__(self, dir_: Path) -> None:
        """
        Initializes this workflow.
        :param dir_: Working directory
        :return: None
        """
        self._dir = dir_
        if not self._dir.exists():
            logging.info(f'Creating working directory: {self._dir}')
            self._dir.mkdir(parents=True)
        self._informs = []

    def run(self, fasta_in: Path, vcf_in: Path, name: str, prefix: str = 'consensus', description: str = '') -> \
            ApplyVariantsOutput:
        """
        Runs the variant caller.
        :param fasta_in: Original FASTA file
        :param vcf_in: Input VCF file
        :param name: Dataset name for sequence
        :param prefix: Prefix for the output file
        :param description: Description for the headers
        :return: None
        """
        path_fasta_out = self._dir / f'{prefix}.fasta'
        with tempfile.TemporaryDirectory(prefix='camel_', dir=Camel.get_instance().config['temp_dir']) as dir_temp:
            path_vcf_idx = self.__index_vcf(vcf_in, Path(dir_temp))
            self.__apply_variants(fasta_in, path_vcf_idx, path_fasta_out)
            self.__update_headers(path_fasta_out, name, description)
        return ApplyVariantsOutput(path_fasta_out, self._informs)

    def __index_vcf(self, path_vcf_in: Path, dir_: Path) -> Path:
        """
        Creates an indexed VCF file.
        :param path_vcf_in: Input VCF file
        :param dir_: Working directory
        :return: Path to indexed VCF file
        """
        path_vcf_out = Path(dir_, f'{path_vcf_in.name}.gz')
        command = Command(' '.join([
            'module load bcftools/1.17;',
            f'bcftools view --output-type z --apply-filters "PASS,." {path_vcf_in} > {path_vcf_out};',
            f'bcftools index -f {path_vcf_out};'
        ]))
        command.run(dir_)
        return path_vcf_out

    def __apply_variants(self, fasta_in: Path, vcf_in: Path, fasta_out: Path) -> None:
        """
        Applies the variants to the target FASTA file.
        :param fasta_in: Input FASTA file
        :param vcf_in: Input VCF file
        :param fasta_out: Output FASTA file
        :return: Path to updated consensus
        """
        command = Command(' '.join([
            'module load bcftools/1.17;'
            f'bcftools consensus --fasta-ref {fasta_in} {vcf_in} > {fasta_out};'
        ]))
        command.run(self._dir)
        if not command.returncode == 0:
            raise RuntimeError(f'Error applying variants: {command.stderr}')
        self._informs.append({'_name_full': 'bcftools consensus 1.17', '_version': '1.17', '_command': command.command})

    def __update_headers(self, path_fasta: Path, name: str, description: str = '') -> None:
        """
        Updates the headers of the FASTA file by adding the target description to the headers.
        :param path_fasta: FASTA file to update
        :param name: Dataset name
        :param description: Description to add
        :return: None
        """
        with open(path_fasta) as handle:
            seqs = list(SeqIO.parse(handle, 'fasta'))
        for s in seqs:
            s.id = f"{name}-{s.id.split('-')[-1]}"
            s.description = description
        with path_fasta.open('w') as handle:
            SeqIO.write(seqs, handle, 'fasta')
