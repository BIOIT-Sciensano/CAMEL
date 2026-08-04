from pathlib import Path
from typing import Any

from camel.snakefiles import human_read_scrubbing

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

INPUT_FASTA_REF = 'iterative_mapping/input/fasta.io'
INPUT_FASTQ = 'iterative_mapping/input/fastq.io'
OUTPUT_REPORT = 'iterative_mapping/report/html.iob'
OUTPUT_SUMMARY = 'iterative_mapping/summary/iterative_mapping.{ext}'
OUTPUT_INFORMS = 'iterative_mapping/report/informs.io'
OUTPUT_FASTA_CONSENSUS_FINAL = 'iterative_mapping/output/fasta.io'
OUTPUT_FASTA_CONSENSUS_FINAL_TRIMMED = 'iterative_mapping/trim_edges/fasta-trim.io'


def get_fasta(config: dict[str, Any]) -> Path:
    """
    Returns the consensus sequence output IO object path.
    """
    if (config['input']['type'] == 'fasta') and ('human_read_scrubbing' not in config['analyses_selected']):
        return Path(str(human_read_scrubbing.INPUT_FASTA).format(input_format='fasta'))
    if (config['input']['type'] == 'fasta') and ('human_read_scrubbing' in config['analyses_selected']):
        return Path(str(human_read_scrubbing.OUTPUT_FASTA).format(input_format='fasta'))
    if config['input']['type'] in ['illumina', 'ont', 'hybrid']:
        return Path(OUTPUT_FASTA_CONSENSUS_FINAL_TRIMMED)
    raise ValueError(f"Invalid input type: {config['input']['type']}")
