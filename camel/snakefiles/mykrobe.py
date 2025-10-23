from pathlib import Path
from typing import Any

from camel.snakefiles import assembly

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

OUTPUT_REPORT = 'mykrobe/report/html.iob'
OUTPUT_REPORT_EMPTY = 'mykrobe/report/html-empty.iob'
OUTPUT_SUMMARY = 'mykrobe/summary_mykrobe.{ext}'
OUTPUT_INFORMS = 'mykrobe/tool/informs.io'


def get_input(config: dict[str, Any]) -> str:
    """
    Returns the input for the Mykrobe tool.
    :param config: Snakemake configuration
    :return: Path to the input file
    """
    if config['input_type'] in ('illumina', 'ont'):
        return 'fq_dict.io'
    if config['input_type'] == 'fasta':
        return assembly.OUTPUT_FASTA
    raise ValueError(f'Unsupported input type: {config["input_type"]}')
