import gzip
import logging
from pathlib import Path
from typing import Optional

import shutil
import socket

from camel.app.config import config
from camel.app.core.utils import fileutils
from camel.app.scriptutils import mainscriptutils


def store_config_file(
        config_file: Path, basename: str, galaxy_job_id: Optional[str] = None, dir_: Path = None) -> Path:
    """
    Exports the config file from the pipeline.
    :param config_file: Config file to export
    :param basename: Basename for tool / pipeline
    :param galaxy_job_id: (Optional) Galaxy job id
    :param dir_: (Optional) Directory to store file, defaults to value from config
    :return: Path to config file
    """
    # Determine the output directory
    if dir_ is not None:
        output_dir = dir_
    else:
        dir_out = config.dir_configs
        if dir_out is None:
            raise RuntimeError("Logging enabled but 'dir_configs' not set in the config.")
        output_dir = Path(dir_out)
    if not output_dir.exists():
        raise FileNotFoundError(f'Logging directory does not exist: {output_dir}')

    # Determine the output path
    export_path = output_dir / '{}.yml.gz'.format('__'.join([
        'config',
        fileutils.make_valid(basename),
        fileutils.make_valid(socket.gethostname()),
        mainscriptutils.get_timestamp_str(),
        galaxy_job_id if galaxy_job_id is not None else 'NA',
    ]))

    # Save the compressed config file
    with gzip.open(export_path, 'wb') as file_out, config_file.open('rb') as file_in:
        shutil.copyfileobj(file_in, file_out)
    logging.info(f"Config file stored in: '{export_path}'")
    return export_path


def store_log_file(log_file: Path, basename: str, galaxy_job_id: Optional[str] = None, is_error_log: bool = False,
                   dir_: Optional[Path] = None) -> Path:
    """
    Stores a log file on disk.
    :param log_file: Log file to store
    :param basename: Basename
    :param galaxy_job_id: (Optional) Galaxy job id
    :param is_error_log: Boolean to indicate if this is an error log
    :param dir_: (Optional) Directory to store file, defaults to value from config
    :return: Path to log file
    """
    # Determine the output directory
    key_dir = 'dir_error_logs' if is_error_log else 'dir_camel_logs'
    # TODO FIX!
    output_dir = dir_ if dir_ is not None else config.dir_logs
    if not output_dir.exists():
        raise RuntimeError(f'Logging directory does not exist: {output_dir}')

    # Determine output file
    prefix = 'error' if is_error_log else 'camel'
    output_path = output_dir / '{}.txt.gz'.format('__'.join([
        prefix,
        fileutils.make_valid(basename).lower(),
        fileutils.make_valid(socket.gethostname()),
        mainscriptutils.get_timestamp_str(),
        galaxy_job_id if galaxy_job_id is not None else 'NA'
    ]))

    # Save log file
    with gzip.open(output_path, 'wb') as file_out, log_file.open('rb') as file_in:
        shutil.copyfileobj(file_in, file_out)
    logging.debug(f"Log file stored in: '{output_path}'")
    return output_path
