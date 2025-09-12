from pathlib import Path
from typing import Any

from camel.resources.snakefile import read_simulation

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_REPORT = 'lrefinder/report/html.iob'
OUTPUT_REPORT_EMPTY = 'lrefinder/report/html_empty.iob'
OUTPUT_INFORMS = 'lrefinder/tool/informs.io'
OUTPUT_SUMMARY = 'lrefinder/summary/summary.{ext}'


def get_input(config: dict[str, Any]) -> str:
    """
    Returns the input for the LREFinder tool.
    :param config: Snakemake configuration
    :return: Path to the input file
    """
    if config['input_type'] in ('illumina', 'ont', 'hybrid'):
        return 'fq_dict.io'
    if config['input_type'] == 'fasta':
        return read_simulation.OUTPUT_FASTQ
    raise ValueError(f"Input type '{config['input_type']}' is not supported by LRE-finder")
