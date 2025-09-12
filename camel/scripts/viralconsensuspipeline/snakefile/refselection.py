import re
from pathlib import Path

SNAKEFILE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

OUTPUT_REPORT = 'ref_selection/report/html.iob'
OUTPUT_REPORT_EMPTY = 'ref_selection/report/html-empty.iob'
OUTPUT_FASTA = 'ref_selection/create_fasta/fasta.io'
OUTPUT_SUMMARY = 'ref_selection/summary/summary.{ext}'


def get_segments(path_db: Path) -> list[str]:
    """
    Returns the names of the segments in the current database.
    :param path_db: Path to the database
    """
    segments = []
    for file_ in (path_db / 'mash').iterdir():
        if not file_.name.endswith('.msh'):
            continue
        segments.append(re.search(r'(.*)\.msh', file_.name).group(1))
    return sorted(segments)
