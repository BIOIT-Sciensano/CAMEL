import dataclasses
from pathlib import Path

from camelcore.app.utils import fastautils, fileutils

from camel.app.loggers import logger
from camel.app.scriptutils import model


@dataclasses.dataclass(frozen=True)
class FastaInput(model.BaseInput):
    """
    Tool FASTA input.
    """
    fasta: Path = dataclasses.field(metadata={'help': 'Input FASTA file'})
    fasta_name: str | None = dataclasses.field(default=None, metadata={'help': 'Name of the input FASTA file'})

    @property
    def name(self) -> str:
        """
        Returns the input name.
        :return: Input name
        """
        name = self.fasta_name if self.fasta_name is not None else self.fasta.name
        for ext in (".fna", ".fasta", ".fa"):
            if name.endswith(ext):
                name = name[:-len(ext)]
                break
        return name

    @property
    def input_str(self) -> str:
        """
        Returns the input file string for the output report.
        :return: Input string
        """
        if self.fasta_name is not None:
            return self.fasta_name
        return self.fasta.name

    def create_symlinks(self, dir_: Path) -> 'FastaInput':
        """
        Creates symlinks for the input files.
        :param dir_: Directory to create symlinks in
        :return: None
        """
        dir_.mkdir(exist_ok=True, parents=True)
        path_symlink = dir_ / (fileutils.make_valid(self.fasta_name) if self.fasta_name is not None else self.fasta.name)
        path_symlink.symlink_to(self.fasta)
        logger.debug(f'Created symlink for FASTA input: {path_symlink} -> {self.fasta}')
        return dataclasses.replace(self, fasta=path_symlink)

    def validate(self) -> bool:
        """
        Checks if the input is valid.
        :return: True if valid, False otherwise
        """
        if fastautils.count_reads(self.fasta) == 0:
            raise ValueError(f'Input FASTA file {self.fasta} is empty.')
        if fastautils.has_duplicates(self.fasta):
            raise ValueError(f'Input FASTA file {self.fasta} contains duplicate sequence IDs.')
        return True
