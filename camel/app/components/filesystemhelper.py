import binascii
import datetime
import re
from pathlib import Path

from camel.app.command.command import Command
from camel.app.loggers import logger


class FileSystemHelper:
    """
    This class contains utility functions to work with the file system.
    """

    TIMESTAMP_FILENAME = "%Y%m%d-%H%M%S"

    @staticmethod
    def make_valid(value: str) -> str:
        """
        Converts arbitrary strings to URL- and filename friendly values.
        :param value: Input value
        :return: URL- and filename friendly value
        """
        value = value.replace(' ', '_')
        return "".join([c for c in value if re.match(r'[\w\-_\\.]', c)])

    @staticmethod
    def get_file_with_extension(input_folder: Path, extension: str) -> Path:
        """
        Returns a single file with the given extension from the given directory.
        :param input_folder: Input directory
        :param extension: File extension
        :return: Path to file
        """
        all_files = FileSystemHelper.get_files_with_extension(input_folder, extension)
        if len(all_files) == 0:
            raise FileNotFoundError(f"No {extension} file found in '{input_folder}'")
        elif len(all_files) > 1:
            raise FileNotFoundError(f"Multiple {extension} files found in '{input_folder}'")
        return all_files[0]

    @staticmethod
    def get_files_with_extension(folder: Path, extension: str) -> list[Path]:
        """
        Returns the files with the given extension from the folder.
        :param folder: Input folder
        :param extension: File extension
        :return: List of paths to files
        """
        return [file_ for file_ in folder.iterdir() if file_.suffix.endswith(extension)]

    @staticmethod
    def get_timestamp_str(timestamp: datetime.datetime = datetime.datetime.now()) -> str:
        """
        Returns the given time stamp as a string that can be used in a filename.
        :param timestamp: Timestamp (default to current time)
        :return: Timestamp as string
        """
        return timestamp.strftime(FileSystemHelper.TIMESTAMP_FILENAME)

    @staticmethod
    def is_gzipped(path: Path) -> bool:
        """
        Checks if the given file is compressed with gzip.
        :param path: Path
        :return: True if gzipped, False otherwise
        """
        with path.open('rb') as handle:
            magic_number = binascii.hexlify(handle.read(2))
        return magic_number == b'1f8b'

    @staticmethod
    def gzip_extract(input_gz_file: Path, output_gz_file: Path) -> None:
        """
        Extracts a GZIP compressed file, the original file is left untouched.
        :param input_gz_file: Input GZ file
        :param output_gz_file: Output path
        :return: None
        """
        logger.info(f"Extracting: {input_gz_file}")
        command = Command(f'gunzip -k -c {input_gz_file} > {output_gz_file}')
        command.run(Path.cwd())
        if not command.returncode == 0:
            raise RuntimeError(f"Cannot extract '{input_gz_file}': {command.stderr}")

    @staticmethod
    def gzip_file(input_file: Path, output_gz_file: Path) -> None:
        """
        Extracts a GZIP compressed file, the original file is left untouched.
        :param input_file: Input non GZ file
        :param output_gz_file: Output path
        :return: None
        """
        logger.info(f"Extracting: {input_file}")
        command = Command(f'gzip -c {input_file} > {output_gz_file}')
        command.run(Path.cwd())
        if not command.returncode == 0:
            raise RuntimeError(f"Cannot gzip '{input_file}': {command.stderr}")

    @staticmethod
    def pigz_extract(input_gz_file: Path, output_gz_file: Path, threads: int) -> None:
        """
        Extracts a GZIP compressed file using pigz, the original file is left untouched.
        :param input_gz_file: Input GZ file
        :param output_gz_file: Output path
        :param threads: Number of threads to use
        :return: None
        """
        logger.info(f"Extracting: {input_gz_file}")
        command = Command(f'pigz -k -p {threads} -dc {input_gz_file} > {output_gz_file}')
        command.run(Path.cwd())
        if not command.returncode == 0:
            raise RuntimeError(f"Cannot extract '{input_gz_file}': {command.stderr}")

    @staticmethod
    def pigz_file(input_file: Path, output_gz_file: Path, threads: int) -> None:
        """
        Compresses a non GZ file using pigz, the original file is left untouched.
        :param input_file: Input non GZ file
        :param output_gz_file: Output path
        :param threads: Number of threads to use
        :return: None
        """
        logger.info(f"Compressing: {input_file}")
        command = Command(f'pigz -c -p {threads} {input_file} > {output_gz_file}')
        command.run(Path.cwd())
        if not command.returncode == 0:
            raise RuntimeError(f"Cannot gzip '{input_file}': {command.stderr}")
