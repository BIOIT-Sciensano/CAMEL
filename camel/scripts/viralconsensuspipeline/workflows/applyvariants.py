import dataclasses
import tempfile
from pathlib import Path

from Bio import SeqIO
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.config import config
from camel.app.loggers import logger
from camel.app.tools.bcftools.bcftoolsconsensus import BcftoolsConsensus
from camel.app.tools.bcftools.bcftoolsindex import BcftoolsIndex
from camel.app.tools.bcftools.bcftoolsview import BcftoolsView


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
            logger.info(f'Creating working directory: {self._dir}')
            self._dir.mkdir(parents=True)

    def run(self, fasta_in: Path, vcf_in: Path, name: str, prefix: str = 'consensus', description: str = '') -> \
            ApplyVariantsOutput:
        """
        Runs the apply variants workflow.
        :param fasta_in: Original FASTA file
        :param vcf_in: Input VCF file
        :param name: Dataset name for sequence headers
        :param prefix: Prefix for the output file
        :param description: Description for the headers
        :return: Apply variants output
        """
        path_fasta_out = self._dir / f'{prefix}.fasta'
        informs = []
        with tempfile.TemporaryDirectory(prefix='camel_', dir=config.dir_temp) as dir_temp:
            path_vcf_idx = self._index_vcf(vcf_in, Path(dir_temp))
            informs.append(self._apply_variants(fasta_in, path_vcf_idx, path_fasta_out))
            self._update_headers(path_fasta_out, name, description)
        return ApplyVariantsOutput(path_fasta_out, informs)

    def _index_vcf(self, path_vcf_in: Path, dir_: Path) -> Path:
        """
        Creates an indexed VCF file.
        :param path_vcf_in: Input VCF file
        :param dir_: Working directory
        :return: Path to indexed VCF file
        """
        path_vcf_out = Path(dir_, f'{path_vcf_in.name}.gz')
        bcftools_view = BcftoolsView()
        bcftools_view.add_input_files({
            'VCF': [ToolIOFile(path_vcf_in)],
        })
        bcftools_view.update_parameters(
            apply_filters='"PASS,."',
            output_filename=str(path_vcf_out),
            output_type='z'
        )
        bcftools_view.run(path_vcf_out.parent)
        bcftools_index = BcftoolsIndex()
        bcftools_index.add_input_files({
            'VCF_GZ': [ToolIOFile(path_vcf_out)]
        })
        bcftools_index.run(path_vcf_out.parent)
        return path_vcf_out

    def _apply_variants(self, fasta_in: Path, vcf_in: Path, fasta_out: Path) -> dict:
        """
        Applies the variants to the target FASTA file.
        :param fasta_in: Input FASTA file
        :param vcf_in: Input VCF file
        :param fasta_out: Output FASTA file
        :return: Tool informs
        """
        bcftools_consensus = BcftoolsConsensus()
        bcftools_consensus.add_input_files({
            'FASTA': [ToolIOFile(fasta_in)],
            'VCF': [ToolIOFile(vcf_in)]
        })
        bcftools_consensus.update_parameters(output_filename=str(fasta_out))
        bcftools_consensus.run()
        return bcftools_consensus.informs

    def _update_headers(self, path_fasta: Path, name: str, description: str = '') -> None:
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
