from pathlib import Path
from typing import Any

from camel.app.toolkits.viralconsensus import dbutils

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

# Per segment
OUTPUT_INFORMS_MASH = 'db_comparison/mash_dist/{segment}/informs.io'
OUTPUT_INFORMS_ALIGN = 'db_comparison/align/{segment}/informs.io'

# Report and summary
OUTPUT_REPORT = 'db_comparison/report/html.iob'
OUTPUT_REPORT_EMPTY = 'db_comparison/report/html-empty.iob'
OUTPUT_SUMMARY = 'db_comparison/summary/summary.{ext}'


def is_enabled(config: dict[str, Any]) -> bool:
    """
    Checks if the database comparison is enabled (selected, reference database available and reads as input).
    :param config: Pipeline configuration
    :return: True if enabled
    """
    return (
        'db_comparison' in config['analyses_selected'] and
        config['ref_selection'].get('db') is not None and
        config['input']['type'] != 'fasta'
    )


def get_db_path(config: dict[str, Any], segment: str, key: str) -> Path:
    """
    Returns the path of a database file used for the database comparison of the given segment.
    :param config: Pipeline configuration
    :param segment: Segment name
    :param key: Database file key ('fasta', 'mash')
    :return: Path to the database file
    """
    return dbutils.get_db_paths(Path(config['ref_selection']['db']), segment)[key]


def get_segments(checkpoints: Any) -> list[str]:
    """
    Returns the segments included in the database comparison (evaluated after splitting the consensus sequence).
    :param checkpoints: Snakemake checkpoints
    :return: Segment names
    """
    dir_consensus = Path(checkpoints.db_comparison_split_consensus.get().output.DIR)
    return sorted(p.stem for p in dir_consensus.glob('*.fasta'))


def get_command_informs(config: dict[str, Any], checkpoints: Any) -> list[str]:
    """
    Returns the informs of the mash dist and minimap2 steps for the command section (first segment only).
    :param config: Pipeline configuration
    :param checkpoints: Snakemake checkpoints
    :return: Informs files
    """
    if not is_enabled(config):
        return []
    segments = get_segments(checkpoints)
    if len(segments) == 0:
        return []
    return [path.format(segment=segments[0]) for path in (OUTPUT_INFORMS_MASH, OUTPUT_INFORMS_ALIGN)]
