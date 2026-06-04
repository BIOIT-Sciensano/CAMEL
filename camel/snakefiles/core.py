from pathlib import Path
from typing import Any

from camel.snakefiles import (
    human_read_scrubbing,
    read_simulation,
    trimming_illumina,
    trimming_ont,
)

SNAKEFILE = f'{Path(__file__).parent / Path(__file__).stem}.smk'
INPUT_FASTA_IO = Path('input', 'fasta.io')
OUTPUT_SUMMARY_INIT = Path('summary', 'summary-init.{ext}')
OUTPUT_SUMMARY = Path('summary', 'output.{ext}')
OUTPUT_HTML_CITATIONS = Path('report', 'html-citations.iob')


def get_fastq_input_downsampling(config: dict[str, Any], read_key: str) -> str:
    """
    Returns the fastq input for the downsampling step.
    :param config: Snakemake configuration
    :param read_key: fastq_se or fastq_pe
    :return: Path to the input file
    """
    if 'human_read_scrubbing' in config.get('analyses_selected', []):
        return human_read_scrubbing.OUTPUT_FASTQ.format(input_format=read_key)
    else:
        return human_read_scrubbing.INPUT_FASTQ.format(input_format=read_key)


def get_fq_input(input_type: str) -> dict[str, str]:
    """
    Returns the FQ input (which is unpacked into the Snakemake input).
    :param input_type: Input type
    :return: Dictionary with strings as keys and paths as values
    """
    if input_type in ('fasta', 'fasta_with_vcf'):
        return {'FASTQ_PE': read_simulation.OUTPUT_FASTQ}
    elif input_type == 'illumina':
        return {'FASTQ_PE': trimming_illumina.OUTPUT_DICT}
    elif input_type == 'ont':
        return {'FASTQ_SE': trimming_ont.OUTPUT_DICT}
    elif input_type == 'hybrid':
        return {
            'FASTQ_PE': trimming_illumina.OUTPUT_DICT,
            'FASTQ_SE': trimming_ont.OUTPUT_DICT
        }
    else:
        raise ValueError(f'Invalid input type: {input_type}')
