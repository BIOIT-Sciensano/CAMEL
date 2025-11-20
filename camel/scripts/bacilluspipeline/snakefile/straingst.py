from pathlib import Path
from typing import Any

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

# Tool outputs
OUTPUT_HDF5 = 'straingst/{read_type}/tool/hdf5-straingst.io'
OUTPUT_STATS = 'straingst/{read_type}/tool/straingst-stats.io'
OUTPUT_STRAINS = 'straingst/{read_type}/tool/straingst-strains.io'

# Report and summary
OUTPUT_INFORMS = 'straingst/{read_type}/tool/informs.io'
OUTPUT_REPORT = 'straingst/{read_type}/report/html.iob'
OUTPUT_REPORT_EMPTY = 'straingst/{read_type}/report/html-empty.iob'
OUTPUT_SUMMARY = 'straingst/{read_type}/summary/summary_out.{ext}'


def get_command_informs(config: dict[str, Any]) -> list[Path]:
    """
    Returns a list of informs IO files for the straingst steps.
    :param config: Snakemake configuration
    :return: List of informs IO files
    """
    paths = []

    # Straingst is disabled -> return empty list
    if not any([an for an in config['analyses'] if an in ['straingst', 'gmo']]):
        return []

    input_type = config['input']['type']
    if input_type in ('illumina', 'hybrid'):
        paths.append(Path(OUTPUT_INFORMS.format(read_type='illumina')))
    if input_type == 'ont':
        paths.append(Path(OUTPUT_INFORMS.format(read_type='ont')))
    return paths


def get_reports(config: dict[str, Any]) -> list[Path]:
    """
    Returns the path to the strainGST report.
    :param config: Snakemake configuration
    :return: List of paths to the straingst reports
    """
    paths = []
    input_type = config['input']['type']

    # FASTA input
    if config['input']['type'] == 'fasta':
        return []

    # FASTQ input
    if input_type in ('illumina', 'hybrid'):
        if any([an for an in config['analyses'] if an in ['straingst', 'gmo']]):
            paths.append(Path(OUTPUT_REPORT.format(read_type='illumina')))
        else:
            paths.append(Path(OUTPUT_REPORT_EMPTY.format(read_type='illumina')))
    if input_type == 'ont':
        if any([an for an in config['analyses'] if an in ['straingst', 'gmo']]):
            paths.append(Path(OUTPUT_REPORT.format(read_type='ont')))
        else:
            paths.append(Path(OUTPUT_REPORT_EMPTY.format(read_type='ont')))

    # Check if at least one path was added
    if len(paths) == 0:
        raise ValueError(f'Invalid input type: {input_type}')

    return paths


def get_summaries(config: dict[str, Any], ext: str) -> list[Path]:
    """
    Returns the paths to the straingst summary file(s).
    :param config: Snakemake configuration
    :param ext: Extension of the summary files
    :return: List of paths to straingst summary files
    """
    input_type = config['input']['type']

    # FASTA input
    if config['input']['type'] == 'fasta':
        return []

    # FASTQ input
    paths = []
    # StrainGST can be part of the GMO assay or can also be run independently, hence the following conditions cover that
    # scenario.
    if (input_type in ('illumina', 'hybrid')) and any([an for an in config['analyses'] if an in ['straingst', 'gmo']]):
        paths.append(Path(OUTPUT_SUMMARY.format(read_type='illumina', ext=ext)))
    if (input_type == 'ont') and any([an for an in config['analyses'] if an in ['straingst', 'gmo']]):
        paths.append(Path(OUTPUT_SUMMARY.format(read_type='ont', ext=ext)))
    return paths
