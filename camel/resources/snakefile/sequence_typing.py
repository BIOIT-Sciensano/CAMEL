from pathlib import Path
from typing import Any, Optional

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

INPUT_FASTA = 'typing/{scheme}/input/fasta.io'
OUTPUT_REPORT = 'typing/{scheme}/report/html.iob'
OUTPUT_REPORT_EMPTY = 'typing/{scheme}/report/html-empty.iob'
OUTPUT_TSV = 'typing/{scheme}/export_tsv/{locus_type}/tsv.io'
OUTPUT_HITS = 'typing/{scheme}/hits/{locus_type}/hits.iob'
OUTPUT_SUMMARY = 'typing/{scheme}/summary/summary_out.{ext}'
OUTPUT_INFORMS = 'typing/{scheme}/commands/informs.io'
OUTPUT_DB_INFORMS = 'typing/{scheme}/scheme_info/informs.iob'
OUTPUT_ALL_MATCHES = 'typing/{scheme}/tsv_profile_matches.io'


def get_detection_method(config_data: dict[str, Any], scheme_key: str, locus_type: str) -> str:
    """
    Returns the detection method for the given scheme.
    :param config_data: Pipeline config
    :param scheme_key: Scheme key
    :param locus_type: Locus type (DNA / peptide)
    :return: Detection method
    """
    if locus_type == 'peptide':
        return 'blast'
    pipeline_setting = config_data['sequence_typing']['options']['method']
    d = config_data['sequence_typing']['dbs'][scheme_key].get('force_method', pipeline_setting)
    return d

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
