from pathlib import Path
from typing import Any

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

# Input
INPUT_FASTQ = 'human_read_scrubbing/{input_format}/input/fastq.io'
INPUT_FASTA = 'human_read_scrubbing/fasta/input/fasta.io' # input_format='fasta'

# Output
OUTPUT_INFORMS = 'human_read_scrubbing/{input_format}/scrubbing/informs.io'
OUTPUT_REPORT = 'human_read_scrubbing/{input_format}/output/html.iob'
OUTPUT_REPORT_EMPTY = 'human_read_scrubbing/{input_format}/output/html-empty.iob'
OUTPUT_SUMMARY = 'human_read_scrubbing/{input_format}/output/summary_out.{ext}'
OUTPUT_FASTQ = 'human_read_scrubbing/{input_format}/compress/scrubbed/fastq_gz.io'
OUTPUT_FASTA = 'human_read_scrubbing/fasta/output/fasta.io'
OUTPUT_FASTQ_REMOVED = 'human_read_scrubbing/{input_format}/compress/removed/fastq_gz.io'
OUTPUT_FASTA_REMOVED = 'human_read_scrubbing/{input_format}/output/fasta_removed.io'


def get_removed(input_format: str) -> str:
    """
    Returns gzipped fastq(s) of the removed reads or the fasta of the removed contigs that serves as input for the
    HRRT reporter.
    :param input_format: input format, either fastq_pe, fastq_se or fasta
    :return: Path to the removed reads/contigs or None
    """
    if input_format == 'fastq_pe':
        return str(OUTPUT_FASTQ_REMOVED).format(input_format='fastq_pe')
    if input_format == 'fastq_se':
        return str(OUTPUT_FASTQ_REMOVED).format(input_format='fastq_se')
    if input_format == 'fasta':
        return str(OUTPUT_FASTA_REMOVED).format(input_format='fasta')
    raise ValueError(f'Invalid input format for human reads scrubbing: {input_format}')


def get_reports(config: dict[str, Any]) -> list[str]:
    """
    Returns the paths to the human read scrubbing reports.
    :param config: Snakemake configuration
    :return: Report path(s)
    """
    input_type = config['input']['type']
    paths = []

    # FASTA input
    if input_type in ('fasta', 'fasta_with_vcf'):
        if 'human_read_scrubbing' in config['analyses']:
            paths.append(str(OUTPUT_REPORT).format(input_format='fasta'))
        else:
            paths.append(str(OUTPUT_REPORT_EMPTY).format(input_format='fasta'))

    # PE reads
    if input_type in ('illumina', 'hybrid'):
        if 'human_read_scrubbing' in config['analyses']:
            paths.append(str(OUTPUT_REPORT).format(input_format='fastq_pe'))
        else:
            paths.append(str(OUTPUT_REPORT_EMPTY).format(input_format='fastq_pe'))

    # SE reads
    if input_type in ('ont', 'hybrid'):
        if 'human_read_scrubbing' in config['analyses']:
            paths.append(str(OUTPUT_REPORT).format(input_format='fastq_se'))
        else:
            paths.append(str(OUTPUT_REPORT_EMPTY).format(input_format='fastq_se'))

    # Check if at least one path was added
    if len(paths) == 0:
        raise ValueError(f'Invalid input type for human reads scrubbing: {input_type}')

    return paths


def get_command_informs(config: dict[str, Any]) -> list[str]:
    """
    Returns a list of informs IO files for the human read scrubbing steps.
    :param config: Snakemake configuration
    :return: List of informs IO files
    """
    input_type = config['input']['type']
    paths = []

    if 'human_read_scrubbing' not in config['analyses']:
        return []

    # FASTA input
    if (input_type in ('fasta', 'fasta_with_vcf')) and ('human_read_scrubbing' in config['analyses']):
        paths.append(str(OUTPUT_INFORMS).format(input_format='fasta'))

    # PE reads
    if (input_type in ('illumina', 'hybrid')) and ('human_read_scrubbing' in config['analyses']):
        paths.append(str(OUTPUT_INFORMS).format(input_format='fastq_pe'))

    # SE reads
    if (input_type in ('ont', 'hybrid')) and ('human_read_scrubbing' in config['analyses']):
        paths.append(str(OUTPUT_INFORMS).format(input_format='fastq_se'))

    # Check if at least one path was added
    if len(paths) == 0:
        raise ValueError(f'Invalid input type for human reads scrubbing: {input_type}')

    return paths


def get_summaries(config: dict[str, Any], ext: str) -> list[Path]:
    """
    Returns the paths to the human read scrubbing summary file(s).
    :param config: Snakemake configuration
    :param ext: Summary format (TSV / JSON)
    :return: Summary file path(s)
    """
    input_type = config['input']['type']
    paths = []

    if 'human_read_scrubbing' not in config['analyses']:
        return []

    # FASTA input
    if (input_type in ('fasta', 'fasta_with_vcf')) and ('human_read_scrubbing' in config['analyses']):
        paths.append(str(OUTPUT_SUMMARY).format(input_format='fasta', ext=ext))

    # PE reads
    if (input_type in ('illumina', 'hybrid')) and ('human_read_scrubbing' in config['analyses']):
        paths.append(str(OUTPUT_SUMMARY).format(input_format='fastq_pe', ext=ext))

    # SE reads
    if (input_type in ('ont', 'hybrid')) and ('human_read_scrubbing' in config['analyses']):
        paths.append(str(OUTPUT_SUMMARY).format(input_format='fastq_se', ext=ext))

    # Check if at least one path was added
    if len(paths) == 0:
        raise ValueError(f'Invalid input type for human reads scrubbing: {input_type}')

    return paths


def get_output_io(config: dict[str, Any]) -> str:
    """
    Returns the paths to the human read scrubbing output io file(s).
    :param config: Snakemake configuration
    :return: Summary file path(s)
    """
    input_type = config['input']['type']

    # FASTA input
    if (input_type in ('fasta', 'fasta_with_vcf')) and ('human_read_scrubbing' in config['analyses']):
        return str(OUTPUT_FASTA).format(input_format='fasta')

    # PE reads
    if (input_type in ('illumina', 'hybrid')) and ('human_read_scrubbing' in config['analyses']):
        return str(OUTPUT_FASTQ).format(input_format='fastq_pe')

    # SE reads
    if (input_type in ('ont', 'hybrid')) and ('human_read_scrubbing' in config['analyses']):
        return str(OUTPUT_FASTQ).format(input_format='fastq_se')
    raise ValueError(f'Invalid input type: {input_type}')
