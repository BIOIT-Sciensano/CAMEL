from pathlib import Path
from typing import Any

from camel.app.snakemake import snakemakeutils

SNAKEFILE_NEXTCLADE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

# Subtype determination
OUTPUT_SUBTYPE_REPORT = 'nextclade/subtype_determination/report/html.iob'
OUTPUT_SUBTYPE_REPORT_EMPTY = 'nextclade/subtype_determination/report/html-empty.iob'

# Output report
OUTPUT_REPORT = 'nextclade/report/html.iob'
OUTPUT_REPORT_EMPTY = 'nextclade/report/html-empty.iob'

# Informs & summary
OUTPUT_SUMMARY = 'nextclade/summary_nextclade.{ext}'
OUTPUT_INFORMS = 'nextclade/informs.io'
OUTPUT_INFORMS_MASH = 'nextclade/subtype_determination/mash/informs.io'

def get_nextclade_db(wildcards, cps, segment: str, config: dict) -> str:
    """
    Returns the path to the Nextclade database.
    :param wildcards: Snakemake wildcards
    :param cps: Snakemake checkpoints
    :param segment: Segment
    :param config: Snakemake config
    :return: Path to the database
    """
    if config['nextclade'].get('dbs') is not None:
        return config['nextclade']['dbs'][segment]

    # Get data from subtype determination
    # noinspection PyUnresolvedReferences
    path_informs = Path(cps.nextclade3_detect_subtype_report.get().output['INFORMS'])
    informs_subtype = snakemakeutils.load_object(path_informs)
    return informs_subtype['nextclade_dbs'][segment]


def get_informs_subtype(wildcards: Any, cps: Any) -> Path:
    """
    Returns the subtype determination informs.
    :param wildcards: Snakemake wildcards
    :param cps: Snakemake checkpoints
    :return: Path to subtype informs
    """
    path_informs = Path(cps.nextclade3_detect_subtype_report.get().output['INFORMS'])
    return path_informs


def get_nextclade_output(wildcards, cps: Any, key: str, config: dict = None) -> list[str]:
    """
    Aggregates the Nextclade output based on the database information.
    :param wildcards: Snakemake wildcards
    :param cps: Snakemake checkpoints
    :param config: Configuration
    :param key: Output key (TSV / INFORMS)
    :return: List of Nextclade outputs
    """
    # Extract segments
    if config['nextclade'].get('dbs') is not None:
        segments = list(config['nextclade']['dbs'].keys())
    else:
        # noinspection PyUnresolvedReferences
        path_informs = Path(cps.nextclade3_detect_subtype_report.get().output['INFORMS'])
        informs_subtype = snakemakeutils.load_object(path_informs)
        segments = list(informs_subtype['nextclade_dbs'].keys())

    # Determine output
    if key == 'TSV':
        base_output = Path('nextclade', '{segment}', 'tsv.io')
    elif key == 'INFORMS':
        base_output = Path('nextclade', '{segment}', 'informs.io')
    else:
        raise ValueError(f'Invalid key: {key}')

    # Return list of outputs
    return [str(base_output).format(segment=segment) for segment in segments]
