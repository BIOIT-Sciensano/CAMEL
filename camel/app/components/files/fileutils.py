import errno
import fileinput
import gzip
import hashlib
import pickle
from pathlib import Path
from typing import Any

from camel.app.components.filesystemhelper import FileSystemHelper


class FileUtils:
    """
    Helper object to perform common operations on the file system.
    """

    @staticmethod
    def hash_file(file_path: Path, block_size: int = 65536) -> str:
        """
        Creates a hash for the file with a default block size of 65536 and the sha256 algorithm.
        :param file_path: File that needs to be hashed
        :param block_size: Block size to be used
        :return: String of the hash with alphanumeric symbols
        """
        if not file_path.is_file():
            raise FileNotFoundError(f"'{file_path}' is not a file")
        hasher = hashlib.sha256()
        with file_path.open('rb') as file_to_hash:
            for chunk in iter(lambda: file_to_hash.read(block_size), b''):
                hasher.update(chunk)
        return hasher.hexdigest()

    @staticmethod
    def get_all_files(directory_path: Path) -> list[Path]:
        """
        Returns all files in a directory recursively.
        """
        files_list = []
        for entry in directory_path.glob('**/*'):
            if entry.is_file():
                files_list.append(entry)
        return files_list

    @staticmethod
    def hash_directory(path: Path) -> str:
        """
        Creates a hash for a folder with a default block size of 65536 and the sha256 algorithm.
        :param path: Directory path
        :return: String of the hash with alphanumeric symbols
        """
        hasher = hashlib.sha256()
        for file_ in sorted(FileUtils.get_all_files(path)):
            hasher.update(FileUtils.hash_file(file_).encode('ascii'))
        return hasher.hexdigest()

    @staticmethod
    def hash_value(value: Any) -> str:
        """
        Creates a hash for a value.
        :param value: Value
        :return: String of the hash with alphanumeric symbols
        """
        hasher = hashlib.sha256()
        hasher.update(pickle.dumps(value))
        return hasher.hexdigest()

    @staticmethod
    def silent_remove(file_path: Path) -> None:
        """
        Silently remove a file, if file does not exist, capture the error
        :param file_path: file to be removed (complete path)
        """
        try:
            file_path.unlink()
        except OSError as e:
            if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
                raise  # re-raise exception if a different error occured

    @staticmethod
    def concatenate_files(output_path: Path, input_files: list[Path]):
        """
        Concatenate the input files specified into one output file. If the input is gzipped,
        the output will also be a gzipped file.
        :param input_files: input files to be concatenated
        :param output_path: Filename of the output
        :return: None
        """
        def get_hook(file):
            if FileSystemHelper.is_gzipped(file):
                return lambda file_name, mode: gzip.open(file_name, mode='rt')
            else:
                return open

        fin = fileinput.input(input_files, openhook=get_hook(input_files[0]))
        output_fn = gzip.open if FileSystemHelper.is_gzipped(input_files[0]) else open
        with output_fn(output_path, 'wt') as fout:
            for line in fin:
                fout.write(line)
        fin.close()
