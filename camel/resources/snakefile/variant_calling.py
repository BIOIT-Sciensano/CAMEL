from pathlib import Path
from typing import Any

from camel.resources.snakefile import read_simulation

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

OUTPUT_REPORT = 'variant_calling/report/html.iob'
OUTPUT_REPORT_EMPTY = 'variant_calling/report/html-empty.iob'
OUTPUT_SUMMARY = 'variant_calling/summary/summary_out.{ext}'
OUTPUT_DEPTH_INFORMS = 'variant_calling/depth/informs.io'
OUTPUT_DEPTH_TSV = 'variant_calling/depth/tsv.io'
OUTPUT_UNFILTERED_VCF = 'variant_calling/unzip_vcf/vcf.io'
OUTPUT_UNFILTERED_VCF_GZ = 'variant_calling/norm/vcf_gz-indexed.io'
OUTPUT_CONSENSUS = 'variant_calling/consensus/fasta.io'
OUTPUT_FILTERED_VCF = 'variant_calling/unzip_vcf_filtered/vcf.io'
OUTPUT_FILTERED_VCF_GZ = 'variant_calling/filter_zscore/vcf_gz-indexed.io'
OUTPUT_MAPPING_RATE_INFORMS = 'variant_calling/rate/informs.io'
OUTPUT_MAPPING_INFORMS = 'variant_calling/read_mapping/informs.io'
OUTPUT_INFORMS_ALL = 'variant_calling/informs_all.io'


def get_mapping_fq_input(config: dict[str, Any]) -> str:
    """
    Returns the fq input file path needed for mapping.
    :param config: Snakemake configuration
    :return: Path to the mapping fq input file
    """
    # ONT and hybrid were added because otherwise some tests of the MockPipeline fail
    if config['input_type'] in ('illumina', 'ont', 'hybrid'):
        return 'fq_dict.io'
    if config['input_type'] in ('fasta', 'fasta_with_vcf'):
        return read_simulation.OUTPUT_FASTQ
    raise ValueError(f'Input type {config["input_type"]} is not supported.')


def get_bam(config: dict[str, Any]) -> str:
    """
    Returns the BAM output IO object path.
    :param config: Snakemake configuration
    :return: Path to the BAM file
    """
    # ONT and hybrid were added because otherwise some tests of the MockPipeline fail
    if config['input_type'] in ('illumina', 'fasta', 'hybrid'):
        return 'variant_calling/read_mapping/illumina/bam.io'  # OUTPUT_BAM_ILLUMINA
    if config['input_type'] == 'ont':
        return 'variant_calling/read_mapping/ont/bam.io'  # OUTPUT_BAM_ONT
    if config['input_type'] == 'fasta_with_vcf':
        return 'variant_calling/dummy_bam/bam.io'
    raise ValueError(f'Input type {config["input_type"]} is not supported.')


def get_vcf(config: dict[str, Any]) -> str:
    """
    Returns the VCF output IO object path (before filtering).
    :param config: Snakemake configuration
    :return: Path to the unfiltered VCF file
    """
    if config['input_type'] in ('fasta', 'illumina', 'ont', 'hybrid'):
        return OUTPUT_UNFILTERED_VCF
    if config['input_type'] == 'fasta_with_vcf':
        return 'input/vcf.io'
    raise ValueError(f'Input type {config["input_type"]} is not supported.')


def get_vcf_gz(config: dict[str, Any]) -> str:
    """
    Returns the VCF GZ output IO object path (before filtering).
    :param config: Snakemake configuration
    :return: Path to the unfiltered gzipped VCF file
    """
    if config['input_type'] in ('fasta', 'illumina', 'ont', 'hybrid'):
        return OUTPUT_UNFILTERED_VCF_GZ
    if config['input_type'] == 'fasta_with_vcf':
        return 'variant_calling/gzip/vcf_gz.io'
    raise ValueError(f'Input type {config["input_type"]} is not supported.')


def get_reports(config: dict[str, Any]) -> str:
    """
    Returns the path to the variant calling report.
    :param config: Snakemake configuration
    :return: Report path
    """
    input_type = config['input_type']
    if input_type in ('illumina', 'fasta', 'ont'):
        return OUTPUT_REPORT
    if input_type == 'fasta_with_vcf':
        return OUTPUT_REPORT_EMPTY
    raise ValueError(f'Input type {config["input_type"]} is not supported.')


def get_mapping_informs(config: dict[str, Any]) -> str:
    """
    Returns the paths to the variant calling mapping informs.
    :param config: config
    :return: Path(s) to the variant calling mapping informs
    """
    input_type = config['input_type']

    if input_type in ('illumina', 'fasta', 'fasta_with_vcf', 'hybrid'):
        # For hybrid, Illumina reads are used for now
        return 'variant_calling/read_mapping/illumina/informs.io'
    elif input_type == 'ont':
        return 'variant_calling/read_mapping/ont/informs.io'
    raise ValueError(f'Input type {config["input_type"]} is not supported.')

def get_summaries(config: dict[str, Any]) -> list[Path]:
    """
    Returns the paths to the variant calling summary file(s).
    :param config: Snakemake configuration
    :return: Summary path(s)
    """
    input_type = config['input_type']
    paths = []
    if input_type in ('illumina', 'fasta', 'ont'):
        paths.append(OUTPUT_SUMMARY)
    return [p for p in paths]


def get_command_informs(config: dict[str, Any]) -> list[Path]:
    """
    Returns the paths to the variant calling informs.
    :param config: config
    :return: Path(s) to the variant calling informs
    """
    input_type = config['input_type']
    paths = []
    if input_type in ('illumina', 'fasta', 'ont'):
        paths.append(OUTPUT_INFORMS_ALL)
    return [p for p in paths]
