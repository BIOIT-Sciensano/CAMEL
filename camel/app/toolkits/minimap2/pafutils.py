import re
from pathlib import Path

import pandas as pd

PAF_COLUMNS = [
    'query',
    'query_length',
    'query_start',
    'query_end',
    'strand',
    'target',
    'target_length',
    'target_start',
    'target_end',
    'nb_matches',
    'alignment_length',
    'mapq',
]
MAX_TAGS = 30
CS_PATTERN = re.compile(r'(=[A-Za-z]+|\*[a-z]{2}|[+\-][a-z]+)')


def read_paf(path_paf: Path, tags: list[str], primary_only: bool = True) -> pd.DataFrame:
    """
    Reads a PAF file, the requested optional tags are added as columns (missing value if not present).
    :param path_paf: PAF file
    :param tags: Optional tags to extract (e.g. ['cs', 'NM'])
    :param primary_only: If True, only the primary alignments are kept (tag 'tp:A:P')
    :return: PAF records
    """
    if path_paf.stat().st_size == 0:
        return pd.DataFrame(columns=[*PAF_COLUMNS, *tags])
    cols_tags = list(range(MAX_TAGS))
    data_paf = pd.read_table(path_paf, header=None, names=[*PAF_COLUMNS, *cols_tags])
    if primary_only:
        data_paf = data_paf[_get_tag(data_paf[cols_tags], 'tp') == 'P']
    return data_paf[PAF_COLUMNS].assign(**{tag: _get_tag(data_paf[cols_tags], tag) for tag in tags})


def _get_tag(data_tags: pd.DataFrame, tag: str) -> pd.Series:
    """
    Returns the value of the given optional tag for each PAF record (e.g. 'cs:Z:=ACGT' -> '=ACGT').
    :param data_tags: Optional tag columns of the PAF records
    :param tag: Tag name
    :return: Tag values
    """
    return data_tags.apply(
        lambda row: next((str(v).split(':', 2)[2] for v in row.dropna() if str(v).startswith(f'{tag}:')), None), axis=1
    )


def count_differences(cs: str) -> dict[str, int | float]:
    """
    Counts the differences between the target and query from the long cs tag (minimap2 --cs=long). Positions with an
    ambiguous base (all non-ACGT bases are reported as 'n' by minimap2) are not compared. Each insertion / deletion is
    counted as a single difference.
    :param cs: Long cs tag
    :return: Identity (%, over the compared positions), differences, mismatches, indels and compared positions
    """
    matches, mismatches, indels = 0, 0, 0
    for op in CS_PATTERN.findall(cs):
        if op[0] == '=':
            matches += len(op[1:]) - op.lower().count('n')
        elif op[0] == '*':
            if 'n' not in op:
                mismatches += 1
        else:
            indels += 1
    compared = matches + mismatches
    return {
        'identity': 100 * matches / compared if compared > 0 else 0.0,
        'differences': mismatches + indels,
        'mismatches': mismatches,
        'indels': indels,
        'compared_positions': compared,
    }
