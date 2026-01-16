from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_INFORMS = 'confindr/tool/informs.io'
OUTPUT_REPORT = 'confindr/report/html.iob'
OUTPUT_REPORT_EMPTY = 'confindr/report/html-empty.iob'
OUTPUT_SUMMARY = 'confindr/summary/summary_confindr.{ext}'


def get_report(config) -> str:
    """
    Returns the path to the ConFindr report io file.
    :param config: Snakemake configuration
    :return: Path to report
    """
    if ('confindr' not in config['analyses']) or (config['input']['type'] not in ['illumina', 'hybrid']):
        return OUTPUT_REPORT_EMPTY
    return OUTPUT_REPORT


def get_command_informs(config) -> list[str]:
    """
    Returns the path to the ConFindr informs io file.
    :param config: Snakemake configuration
    :return: Path to informs IO object
    """
    if ('confindr' not in config['analyses']) or (config['input']['type'] not in ['illumina', 'hybrid']):
        return []
    return [OUTPUT_INFORMS]


def get_summary(config) -> list[str]:
    """
    Returns the path to the ConFindr summay file.
    :param config: Snakemake configuration
    :return: Path to summary TSV
    """
    if ('confindr' not in config['analyses']) or (config['input']['type'] not in ['illumina', 'hybrid']):
        return []
    return [OUTPUT_SUMMARY]
