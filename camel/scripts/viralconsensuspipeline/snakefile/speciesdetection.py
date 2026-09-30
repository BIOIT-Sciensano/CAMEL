from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

# Report and summary
OUTPUT_REPORT = 'species_detection/report/html.iob'
OUTPUT_SUMMARY = 'species_detection/summary/summary.{ext}'
OUTPUT_INFORMS = 'species_detection/informs/informs.io'


def is_enabled(config: dict) -> bool:
    """
    Returns True if the species was detected before running the pipeline.
    :param config: Snakemake config
    :return: True if enabled
    """
    return config.get('species_detection') is not None
