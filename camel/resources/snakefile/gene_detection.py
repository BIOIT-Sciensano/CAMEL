from pathlib import Path
from typing import Any, Optional

from camel.resources.snakefile import read_simulation

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
SNAKEFILE_BLAST = Path(__file__).parent / f'{Path(__file__).stem}_blast.smk'
SNAKEFILE_KMA = Path(__file__).parent / f'{Path(__file__).stem}_kma.smk'

# Input files and database
GENE_DETECTION_FASTA = 'gene_detection/{db}/db_manager/fasta.io'
GENE_DETECTION_FASTA_CLUSTERED = 'gene_detection/{db}/db_manager/fasta-clust.io'
OUTPUT_DB_INFORMS = 'gene_detection/{db}/db_manager/informs.iob'
INPUT_FASTA = 'gene_detection/{db}/input/fasta.io'

# Generic output paths with a wildcard for the detection method
OUTPUT_HITS_METHOD = 'gene_detection/{db}/{method}/hits.iob'
OUTPUT_INFORMS_METHOD = 'gene_detection/{db}/{method}/informs.io'

# Selected hits and informs for the given database
OUTPUT_ALL_HITS = 'gene_detection/{db}/hit_selection/selected-hits.iob'
OUTPUT_INFORMS = 'gene_detection/{db}/hit_selection/informs.io'
OUTPUT_COLUMNS = 'gene_detection/{db}/report/informs-columns.io'

OUTPUT_TSV_BLAST = 'gene_detection/{db}/hit_filtering/tsv-filtered.io'

# Report and summary outputs
OUTPUT_REPORT = 'gene_detection/{db}/report/html.iob'
OUTPUT_REPORT_EMPTY = 'gene_detection/{db}/report/html-empty.iob'
OUTPUT_SUMMARY = 'gene_detection/{db}/report/summary_out.{ext}'


def get_gene_detection_report(db_key: str, config: dict[str, Any], analysis_name: Optional[str] = None) -> str:
    """
    Returns the report input for the given database key.
    :param db_key: Database key
    :param config: Pipeline config
    :param analysis_name: Analysis name that is checked
    :return: Report input path
    """
    search_key = analysis_name if analysis_name is not None else db_key
    if search_key not in config['analyses']:
        return str(OUTPUT_REPORT_EMPTY).format(db=db_key)
    return str(OUTPUT_REPORT).format(db=db_key)
