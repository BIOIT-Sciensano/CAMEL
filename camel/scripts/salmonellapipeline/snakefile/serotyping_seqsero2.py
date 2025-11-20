import re
from pathlib import Path
from typing import Any


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

    if 'serotype' not in config['analyses']:
        return []

    # FASTA or hybrid input
    if input_type in ('fasta', 'hybrid'):
        paths.append(str(OUTPUT_KMER_INFORMS))

    # PE reads
    if input_type == 'illumina':
        paths.append(str(OUTPUT_KMER_INFORMS))
        paths.append(str(OUTPUT_ALLELE_INFORMS))
        paths.append(str(OUTPUT_KMERREAD_INFORMS))

    # SE reads
    if input_type == 'ont':
        paths.append(str(OUTPUT_KMER_INFORMS))
        paths.append(str(OUTPUT_KMERREAD_INFORMS))
    return paths


def seqsero2_output_parser(seqsero2_file: Path, seqsero2_mode: str) -> list[str]:
    """
    Parses the output file of a SeqSero2 run, uniform over all three modes.
    :param seqsero2_file: path of the output file
    :param seqsero2_mode: mode of the SeqSero2 run, either seqsero2_kmer, seqsero2_allele, or seqsero2_kmerread
    2. List of result string to be written to tsv file
    """
    with seqsero2_file.open('r') as handle:
        tsv_results = handle.readlines()[2:8]
    tsv_results = [re.sub(r'([^ ]) ([^ ])', r'\1_\2', res).strip("\n") for res in tsv_results]
    tsv_results = [res.replace(':\t', '\t') for res in tsv_results]
    tsv_results = [seqsero2_mode + "_" + x for x in tsv_results]
    return tsv_results
