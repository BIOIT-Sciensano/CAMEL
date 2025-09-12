from pathlib import Path
from typing import Any, Optional

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

INPUT_FASTA = 'typing/{scheme}/input/fasta.io'
OUTPUT_REPORT = 'typing/{scheme}/html.iob'
OUTPUT_REPORT_EMPTY = 'typing/{scheme}/html-empty.iob'
OUTPUT_TSV = 'typing/{scheme}/{locus_type}/tabular/tsv.io'
OUTPUT_HITS = 'typing/{scheme}/{locus_type}/{detection_method}/hits.iob'
OUTPUT_SUMMARY = 'typing/{scheme}/summary_out.{ext}'
OUTPUT_INFORMS = 'typing/{scheme}/informs.io'
OUTPUT_ALL_MATCHES = 'typing/{scheme}/tsv_profile_matches.io'


def get_sequence_typing_report(scheme_key: str, config: dict[str, Any], analysis_name: Optional[str] = None) -> str:
    """
    Returns the report input for the given database key.
    :param scheme_key: Database key
    :param config: Pipeline config
    :param analysis_name: Analysis name that is checked
    :return: Report input path
    """
    search_key = analysis_name if analysis_name is not None else scheme_key
    if search_key not in config['analyses']:
        return OUTPUT_REPORT_EMPTY.format(scheme=scheme_key)
    return OUTPUT_REPORT.format(scheme=scheme_key)
