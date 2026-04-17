from pathlib import Path
from typing import Any

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

# Note that the abritAMR command is called 'report'
OUTPUT_RUN_INFORMS = 'abritamr/run/informs.io'
OUTPUT_REPORT_REPORT_INFORMS = 'abritamr/report/informs.io'

# Workflow outputs
OUTPUT_REPORT = 'abritamr/output_report/html.iob'
OUTPUT_REPORT_EMPTY = 'abritamr/output_report/html-empty.iob'
OUTPUT_SUMMARY = 'abritamr/summary/summary_out.{ext}'


def get_command_informs(config: dict[str, Any]) -> list[str]:
    """
    Returns the paths to the AbriTAMR informs.
    :param config: config
    :return: Path(s) to the AbriTAMR informs
    """
    paths = []

    # AbriTAMR is disabled -> return empty list
    if 'abritamr' not in config['analyses_selected']:
        return []
    else:
        paths.append(OUTPUT_RUN_INFORMS)
        paths.append(OUTPUT_REPORT_REPORT_INFORMS)

    return paths
