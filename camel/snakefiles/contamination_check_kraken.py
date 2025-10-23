from pathlib import Path
from typing import Any

from camel.snakefiles import assembly

SNAKEFILE = f'{Path(__file__).parent / Path(__file__).stem}.smk'
OUTPUT_REPORT = 'contamination_check/{input_format}/report/html.iob'
OUTPUT_REPORT_EMPTY = 'contamination_check/{input_format}/report/html-empty.iob'
OUTPUT_INFORMS = 'contamination_check/{input_format}/kraken2/informs-contamination.io'
OUTPUT_INFORMS_KRAKEN2 = 'contamination_check/{input_format}/kraken2/informs.io'
OUTPUT_SUMMARY = 'contamination_check/{input_format}/summary/summary_out.{ext}'


def get_input(config: dict[str, Any]) -> str:
    """
    Returns the input for the contamination check.
    :param config: Snakemake configuration
    :return: Path to the input file
    """
    if config['input_type'] in ('illumina', 'ont', 'hybrid'):
        return 'fq_dict.io'
    if config['input_type'] in ('fasta', 'fasta_with_vcf'):
        return assembly.OUTPUT_FASTA
    raise ValueError(f'Invalid input type: {config["input_type"]}')


def get_reports(config: dict[str, Any]) -> list[str]:
    """
    Returns the paths to the contamination check reports.
    :param config: Snakemake configuration
    :return: Report path(s)
    """
    input_type = config['input_type']
    paths = []

    # PE reads
    if input_type in ('illumina', 'hybrid'):
        if 'kraken2' in config['analyses']:
            paths.append(str(OUTPUT_REPORT).format(input_format='fastq_pe'))
        else:
            paths.append(str(OUTPUT_REPORT_EMPTY).format(input_format='fastq_pe'))

    # SE reads
    if input_type in ('ont', 'hybrid'):
        if 'kraken2' in config['analyses']:
            paths.append(str(OUTPUT_REPORT).format(input_format='fastq_se'))
        else:
            paths.append(str(OUTPUT_REPORT_EMPTY).format(input_format='fastq_se'))

    # FASTA input
    if input_type in ('fasta', 'fasta_with_vcf'):
        if 'kraken2' in config['analyses']:
            paths.append(str(OUTPUT_REPORT).format(input_format='fasta'))
        else:
            paths.append(str(OUTPUT_REPORT_EMPTY).format(input_format='fasta'))
    return paths


def get_summaries(config: dict[str, Any], ext: str) -> list[Path]:
    """
    Returns the paths to the contamination check summary file(s).
    :param config: Snakemake configuration
    :param ext: Summary format (TSV / JSON)
    :return: Report path(s)
    """
    input_type = config['input_type']
    paths = []

    # PE reads
    if (input_type in ('illumina', 'hybrid')) and ('kraken2' in config['analyses']):
        paths.append(str(OUTPUT_SUMMARY).format(input_format='fastq_pe', ext=ext))

    # SE reads
    if (input_type in ('ont', 'hybrid')) and ('kraken2' in config['analyses']):
        paths.append(str(OUTPUT_SUMMARY).format(input_format='fastq_se', ext=ext))

    # FASTA input
    if (input_type in ('fasta', 'fasta_with_vcf')) and ('kraken2' in config['analyses']):
        paths.append(str(OUTPUT_SUMMARY).format(input_format='fasta', ext=ext))
    return paths


def get_command_informs(config: dict[str, Any]) -> list[str]:
    """
    Returns the paths to the Kraken informs.
    :param config: config
    :return: Path(s) to the Kraken informs
    """
    input_type = config['input_type']
    paths = []

    # Kraken 2 is disabled -> return empty list
    if 'kraken2' not in config['analyses']:
        return []

    # PE reads
    if input_type in ('illumina', 'hybrid'):
        paths.append(str(OUTPUT_INFORMS_KRAKEN2).format(input_format='fastq_pe'))

    # SE reads
    if input_type in ('ont', 'hybrid'):
        paths.append(str(OUTPUT_INFORMS_KRAKEN2).format(input_format='fastq_se'))

    if input_type in ('fasta', 'fasta_with_vcf'):
        paths.append(str(OUTPUT_INFORMS_KRAKEN2).format(input_format='fasta'))
    return paths
