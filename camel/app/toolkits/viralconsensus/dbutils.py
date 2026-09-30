import json
from collections import defaultdict
from pathlib import Path

import pandas as pd

from camel.app.loggers import logger


def get_db_paths(dir_db: Path, segment: str) -> dict[str, Path | bool | None]:
    """
    Resolves database file paths for a given segment.

    Checks whether clustered sequences are present in the 'full' subdirectory.
    If present, sequence files are retrieved from 'full'; otherwise, top-level files are used.
    :param dir_db: Database directory
    :param segment: Segment name
    :return: Dictionary with keys 'fasta', 'mash', 'clusters' (None if not available) and 'is_full'
    """
    dir_full = dir_db / 'full'
    is_full = (dir_full / 'fasta_by_segment' / f'{segment}.fasta').exists()
    dir_seqs = dir_full if is_full else dir_db
    path_clusters = dir_db / 'clusters' / f'{segment}.tsv'
    return {
        'fasta': dir_seqs / 'fasta_by_segment' / f'{segment}.fasta',
        'mash': dir_seqs / 'mash' / f'{segment}.msh',
        'clusters': path_clusters if path_clusters.exists() else None,
        'is_full': is_full
    }


def get_segments(dir_db: Path) -> list[str]:
    """
    Returns the names of the segments in the database (based on the mash sketches).
    :param dir_db: Database directory
    :return: Segment names (sorted)
    """
    return sorted(p.stem for p in (dir_db / 'mash').glob('*.msh'))


def read_genome_info(dir_db: Path) -> dict:
    """
    Reads the genome information of the database ('genome_size', optionally 'mutually_exclusive_segments').
    :param dir_db: Database directory
    :return: Genome information
    """
    with (dir_db / 'genome_info.json').open() as handle:
        return json.load(handle)


def read_manifest(dir_db: Path) -> dict:
    """
    Reads the database manifest (database name, download date, ...). Not available for databases created before v3.
    :param dir_db: Database directory
    :return: Manifest data, empty if not available
    """
    path_manifest = dir_db / 'manifest.json'
    if not path_manifest.exists():
        return {}
    with path_manifest.open() as handle:
        return json.load(handle)


def read_metadata(dir_db: Path, accessions: set[str] | None = None) -> pd.DataFrame:
    """
    Reads the sequence metadata of the database. Empty values and '-' are parsed as missing values, all values are
    read as text.
    :param dir_db: Database directory
    :param accessions: If specified, only the metadata of these accessions is returned
    :return: Metadata (accession in the 'id' column)
    """
    path_meta = dir_db / 'sequence_metadata.tsv'
    logger.info(f'Reading sequence metadata from: {path_meta}')
    data_meta = pd.read_table(path_meta, dtype=str, keep_default_na=False, na_values=['-', ''])
    if accessions is not None:
        data_meta = data_meta[data_meta['id'].isin(accessions)]
    return data_meta


def read_seq_lengths(dir_db: Path) -> dict[str, dict[str, int]]:
    """
    Reads the length of the sequences used for the reference selection (from the FASTA index files).
    :param dir_db: Database directory
    :return: Sequence length by sequence id, by segment
    """
    len_by_segment = defaultdict(dict)
    for path_fai in (dir_db / 'fasta_by_segment').glob('*.fai'):
        data_length = pd.read_table(path_fai, usecols=[0, 1], names=['seq_id', 'length'])
        for seq_id, length in data_length.itertuples(index=False, name=None):
            len_by_segment[path_fai.name.split('.')[0]][seq_id] = length
    return len_by_segment
