from pathlib import Path
from typing import Any

from camel.resources.snakefile import read_simulation

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

OUTPUT_REPORT = 'shigatyper/report/html.iob'
OUTPUT_REPORT_EMPTY = 'shigatyper/report/html-empty.iob'
OUTPUT_SUMMARY = 'shigatyper/summary_shigatyper.{ext}'
OUTPUT_INFORMS = 'shigatyper/tool/informs.io'


def get_input(config: dict[str, Any]) -> str:
    """
    Returns the input for the ShigaTyper tool.
    :param config: Snakemake configuration
    :return: Path to the input file
    """
    if config['input_type'] in ('illumina', 'ont'):
        return 'fq_dict.io'
    if config['input_type'] == 'fasta':
        return read_simulation.OUTPUT_FASTQ
    raise ValueError(f"Invalid input type: {config['input_type']}")
