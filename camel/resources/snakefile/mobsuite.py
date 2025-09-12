from pathlib import Path
from typing import Any

from camel.app.loggers import logger
from camel.app.snakemake import snakemakeutils

DIR_MOB_SUITE = Path('mob_suite')
SNAKEFILE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

# Input
INPUT_FASTA = 'mob_suite/input/fasta.io'

# Report and summary
OUTPUT_INFORMS = 'mob_suite/tool/informs.io'
OUTPUT_REPORT = 'mob_suite/report/html.iob'
OUTPUT_REPORT_EMPTY = 'mob_suite/report/html-empty.iob'
OUTPUT_SUMMARY = 'mob_suite/summary/summary_mob_suite.{ext}'

OUTPUT_CONTEXT_REPORT = 'mob_suite/genomic_context/html.iob'
OUTPUT_CONTEXT_REPORT_EMPTY = 'mob_suite/genomic_context/html-empty.iob'

GENOMIC_CONTEXT_DB = {
    'amrfinder': {'key': 'amrfinder', 'title': 'AMRFinder', 'contig': 'Contig id', 'gene': 'Element symbol'},
    'bacmet': {'key': 'bacmet', 'title': 'BacMet', 'contig': 'qseqid', 'gene': 'Gene_name'},
    'gene_detection': {'contig': 'Sequence (read or contig)', 'gene': 'Locus'},
    'resfinder4': {'key': 'resfinder4', 'title': 'ResFinder4', 'contig': 'qseqid', 'gene': 'Gene_name'},
}


def collect_genomic_context_input(snake_in: Any, output_tsv: Path, output_informs: Path) -> None:
    """
    Collects the input for the genomic context tool.
    - For gene detection databases, the title is retrieved from the metadata informs
    - For other supported databases, the title is fixed
    :param snake_in: Snakemake input dictionary
    :param output_tsv: Output TSV IO file
    :param output_informs: Output informs IO file
    :return: None
    """
    input_as_dict = {k: Path(v) if len(v) > 0 else None for k, v in snake_in.items()}

    # Create output
    data_out = []
    for k, io_tsv in input_as_dict.items():
        if k.startswith('INFORMS_') or io_tsv is None:
            continue
        db_type = k.split('_')[1]
        # Gene detection databases
        if db_type == 'gd':
            try:
                # Retrieve title from informs
                path_informs = input_as_dict[k.replace('TSV_', 'INFORMS_')]
                informs = snakemakeutils.load_object(path_informs)
            except KeyError as err:
                logger.error(f"INFORMS for gene detection DB '{k}' are missing")
                raise err
            tsv_out = snakemakeutils.load_object(io_tsv)
            data_out.append({
                'tsv': tsv_out[0] if len(tsv_out) > 0 else None,
                'meta': {
                    'title': informs['title'],
                    'key': k.split('_')[2],
                    **GENOMIC_CONTEXT_DB['gene_detection']}
            })

        # Other databases
        else:
            try:
                data_out.append({
                    'tsv': snakemakeutils.load_object(io_tsv)[0],
                    'meta': GENOMIC_CONTEXT_DB[db_type]
                })
            except KeyError as err:
                logger.error(f"Unknown DB for genomic context: {db_type}")
                raise err

    # Save output
    tsv_out_all = [row['tsv'] for row in data_out if row['tsv'] is not None]
    snakemakeutils.dump_object(tsv_out_all, Path(output_tsv))
    snakemakeutils.dump_object([{
        **row['meta'],
        'idx': tsv_out_all.index(row['tsv']) if row['tsv'] is not None else None
    } for row in data_out], Path(output_informs))
