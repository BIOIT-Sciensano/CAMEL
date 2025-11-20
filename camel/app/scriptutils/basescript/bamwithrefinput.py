import dataclasses
from pathlib import Path

from camel.app.core.utils import sambamutils
from camel.app.loggers import logger
from camel.app.scriptutils import model


@dataclasses.dataclass(frozen=True)
class BAMWithRefInput(model.BaseInput):
    """
    Defines the script input.
    """
    # Input files
    bam: Path = dataclasses.field(metadata={'help': 'Input BAM file'})
    reference: Path = dataclasses.field(metadata={'help': 'Reference genome'})

    # Input file names
    bam_name: str | None = dataclasses.field(default=None, metadata={'help': 'Name of the input BAM file'})
    reference_name: str | None = dataclasses.field(default=None, metadata={'help': 'Reference genome name'})

    @property
    def name(self) -> str:
        """
        Returns the input name.
        :return: Input name
        """
        name = self.bam_name if self.bam_name is not None else self.bam.name
        for ext in (".bam", ".sam"):
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
        return self.bam_name if self.bam_name is not None else self.bam.name

    def create_symlinks(self, dir_: Path) -> 'BAMWithRefInput':
        """
        Creates symlinks for the input files.
        :param dir_: Directory to create symlinks in
        :return: None
        """
        dir_.mkdir(exist_ok=True, parents=True)
        path_symlink = dir_ / (self.bam_name if self.bam_name is not None else self.bam.name)
        path_symlink.symlink_to(self.bam)
        return dataclasses.replace(self, bam=path_symlink)

    def validate(self) -> bool:
        """
        Checks if the input is valid.
        :return: True if valid, False otherwise
        """
        logger.debug('Validating BAM input')
        if sambamutils.get_record_count(self.bam) == 0:
            raise ValueError(f'Input BAM file {self.bam} is empty.')
        return True
