import datetime
import json
import shutil
import tempfile
from pathlib import Path

import yaml
from camelcore.app.utils import fileutils
from pydantic import BaseModel, ValidationError

from camel.app.config import config
from camel.app.loggers import logger
from camel.app.tools.mist.mistdownload import MiSTDownload
from camel.app.tools.mist.mistindex import MiSTIndex


class LocusInfo(BaseModel):
    """
    Locus in a typing scheme.
    """
    allele_id_regex: str | None
    allele_page_url: str | None = None
    fasta_path: str
    fasta_md5: str | None = None
    name: str
    name_sanitized: str
    name_valid: str
    type: str
    url: str | None

class SchemeMetadata(BaseModel):
    """
    Metadata for a typing scheme.
    """
    name: str
    title: str
    origin: str | None = None
    url: str | None = None
    loci: list[LocusInfo]

class DownloadConfig(BaseModel):
    """
    Configuration to download a typing scheme.
    """
    profiles: bool = False
    url: str
    downloader: str


def create_db_upd_file(title: str, download_opts: DownloadConfig, dir_out: Path) -> None:
    """
    Creates the JSON file with DB update information.
    :param title: Database title
    :param download_opts: Download options
    :param dir_out: Output directory
    :return: None
    """
    path_out = dir_out / 'db_update_info.json'
    now = datetime.datetime.now()
    with path_out.open('w') as handle:
        json.dump({
            'title': title,
            'url': download_opts.url,
            'last_update_full': now.isoformat(),
            'last_update_date': now.date().isoformat(),
        }, handle, indent=4, sort_keys=True)

def _mist_download_and_index(dir_: Path, download_opts: DownloadConfig, threads: int) -> tuple[Path, Path | None]:
    """
    Download and indexes a typing database using MiST.
    :param dir_: Output directory
    :param download_opts: Download options
    :param threads: Number of threads (for indexing)
    :return: Path to the indexed database and path to the profiles TSV file (if available)
    """
    # Download
    dir_download = Path(dir_, 'download')
    dir_download.mkdir(exist_ok=True, parents=True)
    mist_download = MiSTDownload()
    mist_download.update_parameters(
        url=str(download_opts.url),
        include_profiles=download_opts.profiles,
        downloader=download_opts.downloader,
        output=str(dir_download))
    mist_download.run(dir_download)

    # Index
    dir_idx = Path(dir_, 'index')
    dir_idx.mkdir(exist_ok=True, parents=True)
    mist_index = MiSTIndex()
    mist_index.update_parameters(output=str(dir_idx), threads=threads)
    mist_index.add_input_files({'TXT': mist_download.tool_outputs['TXT']})
    if download_opts.profiles:
        mist_index.add_input_files({'TSV': mist_download.tool_outputs['TSV']})
    mist_index.run(dir_idx)
    return (
        mist_index.tool_outputs['DIR'][0].path,
        mist_download.tool_outputs['TSV'][0].path if download_opts.profiles else None
    )

def download_db_with_mist(dir_out: Path, metadata: SchemeMetadata, download_opts: DownloadConfig, threads: int) -> None:
    """
    Downloads a sequence typing database using MiST.
    :param dir_out: Output directory
    :param metadata: scheme metadata
    :param download_opts: Download options
    :param threads: Number of threads
    :return: None
    """
    with tempfile.TemporaryDirectory(dir=config.dir_temp, prefix='camel_download_') as dir_temp:
        # Download and index the database
        dir_mist_idx, path_tsv_profiles = _mist_download_and_index(Path(dir_temp), download_opts, threads)

        # Copy the output files
        fileutils.move_directory_contents(dir_mist_idx, dir_out / 'mist')
        if download_opts.profiles:
            shutil.copyfile(path_tsv_profiles, dir_out / 'profiles.tsv')

        # Update the MD5 checksums
        loci_updated = []
        for locus in metadata.loci:
            loci_updated.append(locus.model_copy(update={
                'fasta_md5': fileutils.hash_file(dir_out / 'mist' / locus.fasta_path),
            }))
        metadata_out = metadata.model_copy(update={'loci': loci_updated})

        # Create the scheme metadata file
        create_db_upd_file(title=metadata.title, download_opts=download_opts, dir_out=dir_out)
        with open(dir_out / 'scheme_metadata.json', 'w') as handle:
            json.dump(metadata_out.model_dump(), handle, indent=4, sort_keys=True)
        logger.info(f"Typing DB created: {dir_out}")

def download_from_manifest(path_yml: Path, dir_out: Path, threads: int) -> None:
    """
    Downloads a typing scheme from a manifest file.
    :param path_yml: Manifest file path
    :param dir_out: Output directory
    :param threads: Number of threads
    :return: None
    """
    # Parse the manifest file
    with path_yml.open() as handle_:
        data_in = yaml.safe_load(handle_)
    try:
        scheme_metadata = SchemeMetadata(**data_in['scheme_metadata'])
        download_opts = DownloadConfig(**data_in['download'])
    except ValidationError as err:
        logger.warning(f'Error parsing manifest ({path_yml}): {err}')
        raise err

    # Download the database
    download_db_with_mist(dir_out, scheme_metadata, download_opts, threads=threads)
