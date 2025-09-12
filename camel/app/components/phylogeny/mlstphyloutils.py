import json
from pathlib import Path
from typing import Optional

import pandas as pd

from camel.app.loggers import logger


def is_perfect(record: pd.Series, detection_method: str) -> bool:
    """
    Determines if the given hit is perfect.
    :param record: Input record
    :param detection_method: Allele detection method.
    :return: True if perfect, False otherwise
    """
    if detection_method == 'blast':
        if str(record['Allele']) in ('-', '?'):
            return False
        try:
            identity = float(record['% Identity'])
        except ValueError:
            return False
        if identity != 100.0:
            return False
        len_hsp, len_locus = record['HSP/Locus length'].split('/')
        if not len_hsp == len_locus:
            return False
        return True
    elif detection_method == 'kma':
        if str(record['Allele']) == '-':
            return False
        if float(record['% Identity']) != 100.0:
            return False
        if float(record['% Coverage']) != 100.0:
            return False
        return True
    elif detection_method == 'rapid':
        if str(record['Allele(s)']) == '-':
            return False
        return True
    else:
        raise ValueError(f"Invalid detection method: {detection_method}")


def parse_tsv_typing(tsv_path: Path, detection_method: str, use_temp: bool = True) -> dict[str, str]:
    """
    Parses a tabular output file for the sequence typing assay.
    :param tsv_path: Typing output file
    :param detection_method: Allele detection method
    :param use_temp: If enabled, temporary allele ids are used
    :return: Parsed alleles (key: locus, value: allele as a string)
    """
    allele_data = pd.read_table(tsv_path)
    allele_data['is_perfect_hit'] = allele_data.apply(lambda x: is_perfect(x, detection_method), axis=1)
    # When the hashed TSV file is provided -> check for allele id length
    if use_temp is True:
        if detection_method != 'rapid':
            allele_data['is_perfect_hit'] = allele_data.apply(
                lambda x: x['is_perfect_hit'] or len(str(x['Allele'])) == 6, axis=1)
    col_allele = 'Allele' if detection_method != 'rapid' else 'Allele(s)'
    return {r['Locus']: r[col_allele] if r['is_perfect_hit'] else '-' for _, r in allele_data.iterrows()}


def parse_tsv_typing_list(tsv_in: list[tuple[Path, str]], detection_method: Optional[str] = 'blast',
                          use_temp: bool = True) -> pd.DataFrame:
    """
    Parses a list of tabular input files.
    :param tsv_in: List of input files + file name
    :param detection_method: Detection method for sequence typing
    :param use_temp: If enabled, temporary allele ids are used
    :return: Dictionary of detected alleles by sample name
    """
    sample_names = []
    allele_data = []
    for tabular_file, file_name in tsv_in:
        logger.debug(f'Parsing file: {tabular_file}')
        try:
            isolate_name = Path(file_name).stem
            logger.debug(f'Sample name: {isolate_name}')
        except IndexError:
            raise ValueError(f'Cannot determine sample name from: {file_name}')
        alleles_parsed = parse_tsv_typing(tabular_file, detection_method, use_temp)
        logger.debug(f"{sum(v != '-' for _, v in alleles_parsed.items()):,}/{len(alleles_parsed):,} perfect hits")
        allele_data.append(alleles_parsed)
        sample_names.append(isolate_name)
    return pd.DataFrame(allele_data, index=sample_names, dtype=str)


def __get_typing_output_dir(dir_report: Path, html_key: str) -> Path:
    """
    Gets the typing output directory.
    :param dir_report: Report directory
    :param html_key: HTML key
    :return: Path to typing output directory
    """
    dirs_typing = [d for d in (dir_report / 'sequence_typing').iterdir() if d.is_dir()]

    # Return the directory if there is only one
    if len(dirs_typing) == 1:
        return dirs_typing[0]

    # Return the directory that matches the input key
    try:
        return next(d for d in iter(dirs_typing) if d.name.lower() == html_key.lower())
    except StopIteration:
        dir_names = [dir_.name for dir_ in dirs_typing]
        raise FileNotFoundError(
            f"Scheme with key '{html_key}' not found in: {dir_report} (available schemes: {', '.join(dir_names)})")


def parse_html_typing_list(dirs_in: list[Path], html_key: str, detection_method: Optional[str] = 'blast',
                           use_temp: bool = True) -> pd.DataFrame:
    """
    Parses a list of HTML output directories.
    :param dirs_in: List of input directories
    :param html_key: Key of the typing scheme
    :param detection_method: Detection method for sequence typing
    :param use_temp: If enabled, temporary allele ids are used
    :return: Dictionary of detected alleles by sample name
    """
    sample_names = []
    allele_data = []
    for dir_typing in dirs_in:
        # Retrieve metadata from JSON file
        path_json_meta = __get_typing_output_dir(dir_typing, html_key) / 'sequence_typing.json'
        with path_json_meta.open() as handle:
            typing_metadata = json.load(handle)

        # Add sample data to allele data
        sample_names.append(typing_metadata['sample'])

        # Look for file with hashes
        tsv_typing = None
        if use_temp is True:
            try:
                tsv_typing = next(path_json_meta.parent.glob('typing-*-hashes.tsv'))
                logger.debug(f'TSV file with hashes found: {tsv_typing.name}')
            except StopIteration:
                pass

        # Hashing disabled or no hash allele TSv file found
        if tsv_typing is None:
            tsv_typing = next(path_json_meta.parent.glob('typing-*.tsv'))

        alleles_parsed = parse_tsv_typing(tsv_typing, detection_method)
        logger.debug(f"{sum(v != '-' for _, v in alleles_parsed.items()):,}/{len(alleles_parsed):,} perfect hits")
        allele_data.append(alleles_parsed)
    return pd.DataFrame(allele_data, index=sample_names, dtype=str)
