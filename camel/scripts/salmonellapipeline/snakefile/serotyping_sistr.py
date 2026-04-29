from pathlib import Path
from typing import Any

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

OUTPUT_JSON = 'serotyping/sistr/tool/json.io'
OUTPUT_INFORMS = 'serotyping/sistr/tool/informs.io'
OUTPUT_REPORT = 'serotyping/sistr/report/html.iob'
OUTPUT_REPORT_EMPTY = 'serotyping/sistr/report/html-empty.iob'
OUTPUT_SUMMARY = 'serotyping/sistr/summary/summary_out_sistr.{ext}'


def format_sistr_hit(prediction: dict[str, Any], locus: str, antigen: str) -> str:
    """
    Parse SISTR output for a specific locus and format it as a CSV string.

    The O antigen is derived from the 'serogroup' field, while H antigens
    (h1, h2) use their respective antigen fields.

    :param prediction: Dictionary containing SISTR results for a specific locus.
    :param locus: Locus name (e.g. 'fliC', 'fljB', 'wzx', or 'wzy').
    :param antigen: Antigen identifier ('h1', 'h2', or 'o').
    :return: A comma-separated string of hit properties, or 'missing' if the locus is absent.
    """
    is_missing = prediction['is_missing']
    if is_missing:
        return 'missing'

    hit_properties = [
        locus,
        prediction[antigen if antigen != 'o' else 'serogroup'].replace(',', ';'),
        format(prediction['top_result']['pident'], '.2f'),
        '/'.join(
            [
                str(prediction['top_result']['length']),
                str(prediction['top_result']['qlen']),
            ]
        ),
        prediction['top_result']['stitle'],
        '...'.join(
            [
                str(prediction['top_result']['sstart']),
                str(prediction['top_result']['send']),
            ]
        ),
    ]
    return ','.join(hit_properties)
