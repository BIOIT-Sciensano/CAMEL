from pathlib import Path
from typing import Any

from camel.resources.snakefile import trimming_illumina, trimming_ont


def get_reports(config: dict[str, Any]) -> list[Path]:
    """
    Returns the path to the read trimming report.
    :param config: Snakemake configuration
    :return: List of paths to the trimming reports
    """
    paths = []
    input_type = config['input_type']

    # FASTA input
    if input_type in ('fasta', 'fasta_with_vcf'):
        return []

    # FASTQ input
    if input_type in ('illumina', 'hybrid'):
        paths.append(trimming_illumina.OUTPUT_REPORT)
    if input_type in ('ont', 'hybrid'):
        paths.append(trimming_ont.OUTPUT_REPORT)

    # Check if at least one path was added
    if len(paths) == 0:
        raise ValueError(f'Invalid input type: {input_type}')

    return paths


def get_summaries(config: dict[str, Any]) -> list[Path]:
    """
    Returns the paths to the read trimming summary file(s).
    :param config: Snakemake configuration
    :return: Path to read trimming summary file
    """
    input_type = config['input_type']

    # FASTA input
    if input_type in ('fasta', 'fasta_with_vcf'):
        return []

    # FASTQ input
    paths = []
    if input_type in ('illumina', 'hybrid'):
        paths.append(trimming_illumina.OUTPUT_SUMMARY)
    if input_type in ('ont', 'hybrid'):
        paths.append(trimming_ont.OUTPUT_SUMMARY)
    return paths


def get_command_informs(config: dict[str, Any]) -> list[Path]:
    """
    Returns a list of informs IO files for the read trimming steps.
    :param config: Snakemake configuration
    :return: List of informs IO files
    """
    paths = []
    input_type = config['input_type']
    if input_type in ('illumina', 'hybrid'):
        paths.append(trimming_illumina.OUTPUT_INFORMS)
    if input_type in ('ont', 'hybrid'):
        paths.append(trimming_ont.OUTPUT_INFORMS)
    return paths
