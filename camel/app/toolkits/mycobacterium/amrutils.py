import csv
from enum import Enum
from pathlib import Path
from typing import Any

import pandas as pd
from camelcore.app.reports.htmltablecell import HtmlTableCell

from camel.app.loggers import logger


class ConfidenceLevel(Enum):
    """
    This class is used to represent a resistance level.
    """
    ASSOC_R = 'Associated with R'
    ASSOC_R_int = 'Associated with R (int.)'
    ASSOC_S = 'Not associated with R'
    ASSOC_S_int = 'Not associated with R (int.)'
    UNKNOWN = 'Uncertain significance'
    NOT_IN_DB = 'Not in db.'


CONFIDENCE_LEVEL_DATA = {
    ConfidenceLevel.ASSOC_R: {'color': 'red', 'priority': 5},
    ConfidenceLevel.ASSOC_R_int: {'color': 'red', 'priority': 4},
    ConfidenceLevel.ASSOC_S_int: {'color': 'green', 'priority': 3},
    ConfidenceLevel.ASSOC_S: {'color': 'green', 'priority': 2},
    ConfidenceLevel.UNKNOWN: {'color': 'grey', 'priority': 1},
    ConfidenceLevel.NOT_IN_DB: {'color': 'grey', 'priority': 0}
}


def shorten_nucleotide_str(nucleotide_str, max_length: int = 5) -> str:
    """
    Shortens a nucleotide string when it is necessary.
    :param nucleotide_str: Nucleotide string
    :param max_length: Maximum length
    :return: Shortened string
    """
    return nucleotide_str if len(nucleotide_str) <= max_length else \
        f'{nucleotide_str[:max_length]}... ({len(nucleotide_str)} bp)'


def parse_pileup(path_pileup: Path) -> dict[int, list[int]]:
    """
    Parses the pileup input.
    :param path_pileup: Path to the pileup file
    :return: Dictionary with ACTG counts per position
    """
    logger.info(f'Parsing pileup file: {path_pileup}')
    data_pileup = pd.read_table(
        path_pileup, names=['chr', 'pos', 'ref', 'count', 'read_bases'],
        usecols=range(5), index_col=False, quoting=csv.QUOTE_NONE)
    data_pileup['read_bases'] = data_pileup['read_bases'].apply(lambda x: x.upper())
    logger.info(f'{len(data_pileup):,} positions parsed from pileup')
    actg_counts = {}
    base_list = ('A', 'C', 'T', 'G')
    for row in data_pileup.to_dict('records'):
        counts_ref_allele = [row['read_bases'].count(ref_symbol) for ref_symbol in ('.', ',')]
        counts = [row['read_bases'].count(nt) for nt in base_list]
        counts[base_list.index(row['ref'])] = sum(counts_ref_allele)
        actg_counts[row['pos']] = counts
    return actg_counts


def combine_associations(amr_associations: list[dict[str, Any]]) -> tuple[str, HtmlTableCell]:
    """
    Combines the AMR associations.
    :param amr_associations: AMR associations
    :return: Reformatted
    """
    if len(amr_associations) == 0:
        return '-', HtmlTableCell(ConfidenceLevel.NOT_IN_DB.value, color='grey')
    elif len(set(a['confidence'] for a in amr_associations)) == 1:
        first_association = next(iter(amr_associations))
        confidence_level = CONFIDENCE_LEVEL_DATA[ConfidenceLevel(first_association['confidence'])]
        return ', '.join(a['antibiotic_short'] for a in amr_associations), HtmlTableCell(
            first_association['confidence'], confidence_level['color'])
    else:
        index_by_antibiotic = {
            d: i for i, d in enumerate(sorted(a['antibiotic_short'] for a in amr_associations), start=1)}
        antibiotic_str = ', '.join([f'{d}<sup>{i}</sup>' for d, i in sorted(index_by_antibiotic.items())])

        indices_by_conf_level = {}
        for association in amr_associations:
            if association['confidence'] not in indices_by_conf_level:
                indices_by_conf_level[association['confidence']] = []
            indices_by_conf_level[association['confidence']].append(
                index_by_antibiotic[association['antibiotic_short']])
        conf_str = ', '.join(
            f"{conf}<sup>{','.join(str(x) for x in indices)}</sup>" for conf, indices in indices_by_conf_level.items())
        conf_color = next(iter(sorted(
            [ConfidenceLevel(a['confidence']) for a in amr_associations],
            key=lambda cl: -CONFIDENCE_LEVEL_DATA[cl]['priority'])))

        return antibiotic_str, HtmlTableCell(conf_str, color=CONFIDENCE_LEVEL_DATA[conf_color]['color'])
