from pathlib import Path


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
