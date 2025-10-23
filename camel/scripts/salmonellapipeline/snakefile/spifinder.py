import json
from pathlib import Path
from typing import Any

import pandas as pd

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

OUTPUT_FASTA_JSON = 'spifinder/fasta/spifinder_output.io'
OUTPUT_FASTA_INFORMS = 'spifinder/fasta/informs.io'
OUTPUT_FASTQ_JSON = 'spifinder/spifinder_fastq/spifinder_output.io'
OUTPUT_FASTQ_INFORMS = 'spifinder/spifinder_fastq/informs.io'
OUTPUT_REPORT = 'spifinder/report/html.iob'
OUTPUT_REPORT_EMPTY = 'spifinder/report/html-empty.iob'
OUTPUT_SUMMARY = 'spifinder/summary/summary_out.{ext}'


def get_command_informs(config: dict[str, Any]) -> list[str]:
    """
    Returns a list of informs IO files for SPIFinder.
    :param config: Snakemake configuration
    :return: List of informs IO files
    """
    input_type = config['input_type']
    paths = []

    if 'spifinder' not in config['analyses']:
        return []

    # FASTA input
    if input_type in ('fasta', 'hybrid'):
        paths.append(OUTPUT_FASTA_INFORMS)

    # PE reads
    if input_type == 'illumina':
        paths.append(OUTPUT_FASTA_INFORMS)
        paths.append(OUTPUT_FASTQ_INFORMS)

    # SE reads
    if input_type == 'ont':
        paths.append(OUTPUT_FASTA_INFORMS)
        paths.append(OUTPUT_FASTQ_INFORMS)
    return paths


def get_summaries(config: dict[str, Any]) -> list[str]:
    """
    Returns the paths to the spifinder summary file(s).
    :param config: Snakemake configuration
    :return: Summary file path(s)
    """
    input_type = config['input_type']
    paths = []

    if 'spifinder' not in config['analyses']:
        return []

    # FASTA input
    if input_type == 'fasta':
        paths.append(OUTPUT_SUMMARY)

    # PE reads
    if input_type in ('illumina', 'hybrid'):
        paths.append(OUTPUT_SUMMARY)

    # SE reads
    if input_type in ('ont', 'hybrid'):
        paths.append(OUTPUT_SUMMARY)

    return paths


def parse_json(path_in: Path, mode: str) -> pd.DataFrame:
    """
    Parses the SPI-finder JSON output.
    :param path_in: Path to the input JSON path
    :param mode: FASTA or FASTQ
    :return: DataFrame with detected hits
    """
    with path_in.open() as handle:
        json_file = json.load(handle)

    # Extract the header:
    header = ['SPI', 'identity', 'coverage']
    if mode == 'fasta':
        header.extend(['contig_name', 'positions_in_contig', 'accession', 'insertion_site', 'category_function'])
    elif mode == 'fastq':
        header.extend(['accession', 'insertion_site', 'category_function'])
    else:
        raise ValueError(f"This function's parameter 'mode' must be either fastq or fasta, current value is {mode}")

    # Create the results
    data_spi_as_dict = json_file['spifinder']['results']['Salmonella Pathogenicity Islands']['SPI']
    if data_spi_as_dict == "No hit found":
        return pd.DataFrame(columns=header)
    return pd.DataFrame([v for _, v in data_spi_as_dict.items()])[header]
