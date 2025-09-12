from pathlib import Path
from typing import Any

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

OUTPUT_JSON = 'serotyping/sistr/tool/json.io'
OUTPUT_INFORMS = 'serotyping/sistr/tool/informs.io'
OUTPUT_REPORT = 'serotyping/sistr/report/html.iob'
OUTPUT_REPORT_EMPTY = 'serotyping/sistr/report/html-empty.iob'
OUTPUT_SUMMARY = 'serotyping/sistr/summary/summary_out_sistr.{ext}'


def sistr_output_parser(prediction: dict[str, Any], locus: str, antigen: str, hits_dict_tsv: dict[str, str]) -> None:
    """
    Parses the SISTR output for a specific locus (the o antigen has 2 loci, and the h antigens each have one locus).
    Updates the hits dictionaries in place without returning any output.
    :param prediction: the dictionary of the results of the specific locus
    :param locus: locus name, either fliC, fljB, wzx, or wzy
    :param antigen: antigen name, either h1, h2, or o
    :param hits_dict_tsv: the dictionary of the results for the tsv file
    :return: None
    """
    is_missing = prediction['is_missing']
    if not is_missing:
        hit_properties = [
            locus,
            prediction[antigen if antigen != 'o' else 'serogroup'].replace(',', ';'),
            format(prediction['top_result']['pident'], '.2f'),
            '/'.join([str(prediction['top_result']['length']),
                      str(prediction['top_result']['qlen'])]),
            prediction['top_result']['stitle'],
            '...'.join([str(prediction['top_result']['sstart']),
                        str(prediction['top_result']['send'])])
        ]

        hits_dict_tsv[f'hits_serotype_{antigen}_{locus}'] = ','.join(hit_properties)
