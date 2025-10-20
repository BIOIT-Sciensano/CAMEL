from pathlib import Path
from typing import Any

from camel.resources.snakefile import human_read_scrubbing

SNAKEFILE_ITERATIVE_MAPPING = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_iterative_mapping = Path('iterative_mapping')
INPUT_FASTA_REF = _dir_iterative_mapping / 'input' / 'fasta.io'
INPUT_FASTQ = _dir_iterative_mapping / 'input' / 'fastq.io'
OUTPUT_ITERATIVE_MAPPING_REPORT = _dir_iterative_mapping / 'report' / 'html.io'
OUTPUT_ITERATIVE_MAPPING_SUMMARY = _dir_iterative_mapping / 'report' / 'summary_iterative_mapping.tsv'
OUTPUT_ITERATIVE_MAPPING_INFORMS = _dir_iterative_mapping / 'report' / 'informs.io'
OUTPUT_ITERATIVE_MAPPING_FASTA_CONSENSUS_FINAL = _dir_iterative_mapping / 'output' / 'fasta.io'
OUTPUT_ITERATIVE_MAPPING_FASTA_CONSENSUS_FINAL_TRIMMED = _dir_iterative_mapping / 'trim_edges' / 'fasta-trim.io'


def get_fasta(config: dict[str, Any]) -> Path:
    """
    Returns the consensus sequence output IO object path.
    """
    if (config['input_type'] == 'fasta') and ('human_read_scrubbing' not in config['analyses']):
        return Path(str(human_read_scrubbing.INPUT_SCRUBBING_FASTA).format(input_format='fasta'))
    if (config['input_type'] == 'fasta') and ('human_read_scrubbing' in config['analyses']):
        return Path(str(human_read_scrubbing.OUTPUT_SCRUBBING_FASTA).format(input_format='fasta'))
    if config['input_type'] in ['illumina', 'ont', 'hybrid']:
        return OUTPUT_ITERATIVE_MAPPING_FASTA_CONSENSUS_FINAL_TRIMMED
    raise ValueError(f"Invalid input type: {config['input_type']}")
