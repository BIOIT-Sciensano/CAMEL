from pathlib import Path
from typing import Any

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

INPUT_FASTQ = 'downsampling/{read_key}/input/fastq.io'
OUTPUT_REPORT = 'downsampling/{read_key}/report/html.iob'
OUTPUT_SUMMARY = 'downsampling/{read_key}/summary/summary_downsampling.{ext}'
OUTPUT_FASTQ = 'downsampling/{read_key}/seqtk/fastq.io'
OUTPUT_INFORMS = 'downsampling/{read_key}/seqtk/informs.io'


def get_reports(config: dict[str, Any], input_type: str = None) -> list[Path]:
    """
    Returns the paths to the downsampling reports.
    :param input_type: input type
    :param config: Snakemake configuration
    :return: Report path(s)
    """
    if input_type is None:
        input_type = config['input_type']

    # FASTA input
    if input_type in ('fasta', 'fasta_with_vcf'):
        return []

    # FASTQ input
    paths = []
    if input_type in ('illumina', 'hybrid'):
        paths.append(str(OUTPUT_REPORT).format(read_key='fastq_pe'))
    if input_type in ('ont', 'hybrid', 'iontorrent'):
        paths.append(str(OUTPUT_REPORT).format(read_key='fastq_se'))

    # Check if at least one path was added
    if len(paths) == 0:
        raise ValueError(f'No downsampling report(s) for input type: {input_type}')

    return paths


def get_command_informs(config: dict[str, Any], input_type: str = None) -> list[Path]:
    """
    Returns a list of informs IO files for the downsampling steps.
    :param config: Snakemake configuration
    :param input_type: Input type
    :return: List of informs IO files
    """
    if input_type is None:
        input_type = config['input_type']

    # FASTA input
    if input_type in ('fasta', 'fasta_with_vcf'):
        return []

    # FASTQ input
    paths = []
    if input_type in ('hybrid', 'illumina'):
        paths.append(str(OUTPUT_INFORMS).format(read_key='fastq_pe'))
    if input_type in ('iontorrent', 'hybrid', 'ont'):
        paths.append(str(OUTPUT_INFORMS).format(read_key='fastq_se'))

    # Check if at least one path was added
    if len(paths) == 0:
        raise ValueError(f'No downsampling informs for input type: {input_type}')

    return paths


def get_summaries(config: dict[str, Any], ext: str, input_type: str = None) -> list[Path]:
    """
    Returns the paths to the downsampling summary file(s).
    :param ext: Summary format (TSV / JSON)
    :param input_type: input type
    :param config: Snakemake configuration
    :return: Summary file path(s)
    """
    if input_type is None:
        input_type = config['input_type']

    # FASTA input
    if input_type in ('fasta', 'fasta_with_vcf'):
        return []

    # FASTQ input
    paths = []
    if input_type in ('hybrid', 'illumina'):
        paths.append(str(OUTPUT_SUMMARY).format(read_key='fastq_pe', ext=ext))
    if input_type in ('iontorrent', 'hybrid', 'ont'):
        paths.append(str(OUTPUT_SUMMARY).format(read_key='fastq_se', ext=ext))

    # Check if at least one path was added
    if len(paths) == 0:
        raise ValueError(f'No downsampling summary output for input type: {input_type}')

    return paths
