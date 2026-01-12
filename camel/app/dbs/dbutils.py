import collections
import tempfile
from enum import Enum
from pathlib import Path
from typing import Any

import humanize
import requests
from furl import furl
from pydantic import BaseModel

from camel.app.config import config
from camel.app.core import cameltesthelper
from camel.app.core.utils import fileutils
from camel.app.dbs import typingdbdownload
from camel.app.loggers import logger


class DBMethod(Enum):
    """
    Methods to download DBs.
    """
    FTP = 'ftp'
    MIST = 'mist'
    SKIP = 'skip'


class DBEntry(BaseModel):
    """
    Represents a database entry.
    """
    location: Path
    is_file: bool = False
    required: bool = False
    method: DBMethod = DBMethod.FTP
    download_opts: dict[str, Any] | None = None


def download_file(url: str, local_path: Path, chunk_size: int = 8192) -> None:
    """
    Downloads a file over HTTPS (or HTTP) using streaming.
    :param url: HTTPS/HTTP URL
    :param local_path: Local output file path
    :param chunk_size: Streaming chunk size in bytes
    :return: None
    """
    local_path = Path(local_path)
    local_path.parent.mkdir(parents=True, exist_ok=True)

    with requests.get(url, stream=True) as r:
        # Check if the URL is valid
        r.raise_for_status()
        with local_path.open('wb') as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)


def download_db_from_archive(path_url: str, dir_out: Path) -> None:
    """
    Downloads a CAMEL database from the given URL and extracts it to the given output directory.
    :param path_url: URL to the database archive
    :param dir_out: Output directory
    :return: None
    """
    with tempfile.TemporaryDirectory(dir=config.dir_temp, prefix='camel_download_') as dir_temp:
        # Download the archive
        path_out = Path(dir_temp, furl(path_url).path.segments[-1])
        path_out.parent.mkdir(exist_ok=True, parents=True)
        try:
            download_file(path_url, path_out)
        except Exception as err:
            logger.warning(f'Failed to download: {path_url}')
            raise err
        logger.info(f"File download to: {path_out} ({humanize.naturalsize(path_out.stat().st_size)})")

        # Extract the database
        dir_extract = Path(dir_temp, 'extract')
        dir_extract.mkdir(exist_ok=True, parents=True)
        fileutils.extract_tgz(path_out, dir_extract)
        logger.info(f"Database extracted to: {dir_extract}")

        # Move the extracted files to the target directory
        fileutils.move_directory_contents(dir_extract, dir_out)

def download_dbs(yml: Path, force: bool, keys: list[str] | None, dir_db: Path, threads: int) -> None:
    """
    Downloads databases from the CAMEL FTP server.
    :param yml: Pipeline YML file
    :param force: Force download, even if the DBs already exist
    :param keys: Keys of the databases to download
    :param dir_db: Databases directory
    :param threads: Number of threads to use for indexing
    :return: None
    """
    # Parse the YAML file
    data_dbs = cameltesthelper.extract_from_yaml(yml, 'dbs', placeholders={'DB_ROOT': str(dir_db)})
    target_keys = data_dbs.keys() if keys is None else keys

    # Iterate over the databases
    nb_by_status = collections.Counter()
    for key in target_keys:
        entry = DBEntry(**data_dbs[key])
        relative_path = Path(entry.location).relative_to(dir_db)

        # Check if the database already exists
        if Path(entry.location).exists() and not force:
            logger.info(f"DB '{key}' already exists at: {entry.location}, skipping (use --force to overwrite)")
            nb_by_status['SKIP'] += 1
            continue

        logger.info(f"Downloading DB '{key}' to: {entry.location} (method: {entry.method.value})")
        try:
            if entry.method == DBMethod.FTP:
                download_db_from_archive(
                    path_url=furl(config.ftp_server).add(path='dbs').add(path=f'{relative_path}.tgz').url,
                    dir_out=entry.location.parent if entry.is_file else entry.location,
                )
            elif entry.method == DBMethod.MIST:
                path_url_manifest = furl(config.ftp_server).add(path='dbs').add(path=f'{relative_path}.yml').url
                path_manifest = Path(config.dir_temp, 'manifest.yml')
                logger.info("Downloaded DB manifest")
                download_file(path_url_manifest, path_manifest)
                typingdbdownload.download_from_manifest(path_manifest, dir_out=entry.location, threads=threads)
            elif entry.method == DBMethod.SKIP:
                logger.info(f"The '{key}' database is not available for download, skipping")
                nb_by_status['SKIP'] += 1
                continue
            else:
                raise ValueError(f"Unknown DB method: {entry.method}")
            nb_by_status['OK'] += 1
        except BaseException as err:
            nb_by_status['ERR'] += 1
            logger.warning(f"Failed to download DB '{key}': {err}")
    logger.info(f"Finished downloading databases (OK/SKIP/ERR: {'/'.join(str(nb_by_status.get(x, 0)) for x in ('OK', 'SKIP', 'ERR'))})")
