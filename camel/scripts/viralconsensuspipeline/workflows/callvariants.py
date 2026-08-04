import dataclasses
import logging
from pathlib import Path
from typing import Any

from camelcore.app.io.tooliofile import ToolIOFile
from cyvcf2 import VCF, Variant

from camel.app.core.piping import pipeutils
from camel.app.tools.bcftools.bcftoolscall import BcftoolsCall
from camel.app.tools.bcftools.bcftoolsmpileup import BcftoolsMpileup
from camel.app.tools.clair3.clair3 import Clair3
from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex


@dataclasses.dataclass
class CallVariantsOutput:
    """
    Holder for the output of the variant calling output.
    """
    path_vcf: Path
    stats: dict[str, Any]
    informs: list[dict]


class CallVariants:
    """
    Wrapper around the variant callers for the viral consensus pipeline.
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

    def run(self, bam_in: Path, fasta_ref: Path, input_type: str, caller: str, params: dict[str, Any],
            threads: int = 4) -> CallVariantsOutput:
        """
        Runs the variant caller.
        :param bam_in: Input BAM file
        :param fasta_ref: Reference FASTA file
        :param input_type: Read type
        :param caller: Variant caller (options: GATK, samtools, clair3)
        :param params: Variant calling parameters
        :param threads: Number of threads
        :return: None
        """
        path_ref = self._prepare_ref_genome(fasta_ref)
        if caller == 'bcftools':
            path_vcf = self.call_variants_bcftools(bam_in, path_ref, input_type)
        elif caller == 'clair3':
            path_vcf = self.call_variants_clair3(bam_in, path_ref, input_type, params['model'], threads)
        else:
            raise ValueError(f"Invalid caller: {caller}")
        stats = self._extract_stats(path_vcf)
        return CallVariantsOutput(path_vcf, stats, self._informs)

    def call_variants_bcftools(self, bam_in: Path, path_fasta: Path, input_type: str, max_depth: int = 8_000) -> Path:
        """
        Calls variants using bcftools.
        :param path_fasta: Input FASTA file
        :param bam_in: Input BAM file
        :param input_type: Read type
        :param max_depth: Maximum depth
        :return: Output VCF file
        """
        # Create working directory
        path_vcf_out = self._dir / 'variants.vcf'

        # Pileup
        bcftools_mpileup = BcftoolsMpileup()
        bcftools_mpileup.add_input_files({
            'BAM': [ToolIOFile(bam_in)],
            'FASTA': [ToolIOFile(path_fasta)],
        })
        bcftools_mpileup.update_parameters(
            output_type='v',
            max_depth=max_depth,
            config='illumina' if input_type == 'illumina' else 'ont',
        )

        # Variant calling
        bcftools_call = BcftoolsCall()
        bcftools_call.update_parameters(
            calling_method='consensus',
            output_type='v',
            output_filename=str(path_vcf_out),
            ploidy='1',
            variants_only=True
        )
        pipeutils.run_as_pipe([bcftools_mpileup, bcftools_call], self._dir)
        self._informs.extend([bcftools_mpileup.informs, bcftools_call.informs])
        return path_vcf_out

    def call_variants_clair3(
            self, bam_in: Path, path_fasta: Path, input_type: str, path_model: Path, threads: int) -> Path:
        """
        Calls variants using Clair3.
        :param bam_in: Input BAM file
        :param path_fasta: Input FASTA file
        :param input_type: Read type
        :param path_model: Path to the Clair3 model
        :param threads: Number of threads
        """
        clair3 = Clair3()
        clair3.add_input_files({'FASTA': [ToolIOFile(path_fasta)], 'BAM': [ToolIOFile(bam_in)]})

        # Check if the model exists
        if not path_model.exists():
            raise FileNotFoundError(f'Clair3 model not found: {path_model}')

        # Options that are always enabled
        platform = 'ilmn' if input_type == 'illumina' else 'ont'
        clair3.update_parameters(
            platform=platform, output_path=str(self._dir), threads=threads, haploid_precise=True, include_ctgs=True,
            no_phasing=True, model_path=str(path_model), chunk_size=5_000)

        # Run Clair3
        clair3.run(self._dir)
        self._informs.append(clair3.informs)
        return clair3.tool_outputs['VCF'][0].path

    def _extract_stats(self, path_vcf: Path) -> dict[str, Any]:
        """
        Extracts variant calling stats by parsing the output VCF file.
        :param path_vcf: Input VCF file
        :return: Variant calling statistics
        """
        with VCF(str(path_vcf)) as vcf_reader:
            variants: list[Variant] = list(vcf_reader)
        return {
            'nb_variants': len(variants),
            'nb_snps': sum(v.is_snp for v in variants),
            'nb_indels': sum(v.is_indel for v in variants)
        }

    def _prepare_ref_genome(self, fasta_in: Path) -> Path:
        """
        Prepares the reference genome by creating a symbolic link in the  current directory and creating a samtools
        index.
        :param fasta_in: Input FASTA file
        :return: Path to indexed FASTA file
        """
        # Create symlink
        path_symlink = self._dir / 'ref' / fasta_in.name
        path_symlink.parent.mkdir(exist_ok=True, parents=True)
        if not path_symlink.exists():
            path_symlink.symlink_to(fasta_in)

        # Create FAI index
        samtools_faidx = SamtoolsFastaIndex()
        samtools_faidx.add_input_files({'FASTA': [ToolIOFile(path_symlink)]})
        samtools_faidx.update_parameters(symlink_input=False)
        samtools_faidx.run(path_symlink.parent)
        return path_symlink
