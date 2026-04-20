from pathlib import Path
from typing import Any

from camel.snakefiles import read_simulation

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_KMER_INFORMS = 'serotyping/seqsero2/kmer/informs.io'
OUTPUT_ALLELE_INFORMS = 'serotyping/seqsero2/allele/informs.io'
OUTPUT_KMERREAD_INFORMS = 'serotyping/seqsero2/kmerread/informs.io'

OUTPUT_REPORT = 'serotyping/seqsero2/report/html.iob'
OUTPUT_REPORT_EMPTY = 'serotyping/seqsero2/report/html-empty.iob'
OUTPUT_SUMMARY = 'serotyping/seqsero2/summary/summary_seqsero2.{ext}'


def get_command_informs(config: dict[str, Any]) -> list[str]:
    """
    Returns a list of informs IO files for serotyping.
    :param config: Snakemake configuration
    :return: List of informs IO files
    """
    input_type = config['input']['type']
    paths = []

    if 'serotype' not in config['analyses_selected']:
        return []

    # FASTA or hybrid input
    if input_type in ('fasta', 'hybrid'):
        paths.append(str(read_simulation.OUTPUT_INFORMS))
        paths.append(str(OUTPUT_KMER_INFORMS))

    # PE reads
    if input_type == 'illumina':
        paths.append(str(OUTPUT_KMER_INFORMS))
        paths.append(str(OUTPUT_ALLELE_INFORMS))
        paths.append(str(OUTPUT_KMERREAD_INFORMS))

    # SE reads
    if input_type == 'ont':
        paths.append(str(read_simulation.OUTPUT_INFORMS))
        paths.append(str(OUTPUT_KMER_INFORMS))
        paths.append(str(OUTPUT_ALLELE_INFORMS))
        paths.append(str(OUTPUT_KMERREAD_INFORMS))
    return paths
