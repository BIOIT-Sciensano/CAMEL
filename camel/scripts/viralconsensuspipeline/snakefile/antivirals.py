from pathlib import Path
from typing import Any

from camel.scripts.viralconsensuspipeline.snakefile import nextclade3

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_REPORT = 'antivirals/report/html.iob'
OUTPUT_REPORT_EMPTY = 'antivirals/report/html-empty.iob'
OUTPUT_SUMMARY = 'antivirals/summary/summary.{ext}'


def get_informs_subtype(wildcards: Any, cps: Any, config: dict) -> Path | list:
    """
    Returns the subtype determination informs when the subtype is not specified in the configuration.
    :param wildcards: Snakemake wildcards
    :param cps: Snakemake checkpoints
    :param config: Snakemake config
    :return: Path to subtype informs, empty list if the subtype is specified in the configuration
    """
    if config['antivirals'].get('subtype') is not None:
        return []
    return nextclade3.get_informs_subtype(wildcards, cps)
