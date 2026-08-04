from pathlib import Path

from camelcore.app.io.tooliodirectory import ToolIODirectory

from camel.app.loggers import logger
from camel.app.tools.star.star import Star


class StarIndex(Star):
    """
    Index reference genome using STAR.
    """

    def __init__(self) -> None:
        """
        Initializes STAR 2.7.11b indexing.
        :return: None
        """
        super().__init__('STAR', '2.7.11b')
        self._required_inputs = ['FASTA']
        self._input_string = "--runMode genomeGenerate"

    def _set_input(self) -> None:
        """
        Sets the input specification and the input string.
        :return: None
        """
        option_fasta = "--genomeFastaFiles"
        for fasta in self._tool_inputs['FASTA']:
            input_fasta = self._symlink_fasta(Path(str(fasta))) if 'symlink_input' in self._parameters else (Path(str(fasta)))
            option_fasta += f" {input_fasta}"

        option_gtf = ""
        if 'GTF' in self._tool_inputs:
            option_gtf = f"--sjdbGGTFfile {self._tool_inputs['GTF'][0].path}"

        self._input_string = ' '.join([self._input_string, option_fasta, option_gtf])

    def _set_output(self) -> None:
        """
        Sets the output specification.
        :return: None
        """
        index_dir = self._folder / "GenomeDir" if "symlink_input" in self._parameters else (
                self._tool_inputs['FASTA'][0].path.parent / "GenomeDir")
        self._tool_outputs['INDEX_DIR'] = [ToolIODirectory(index_dir)]

    def _symlink_fasta(self, fasta: Path) -> Path:
        """
        Creates a symlink for the fasta input. This avoids errors when there are no writing permissions
        on the directory of the input fasta.
        :param fasta: Input fasta
        :return: Path to symlink input
        """
        path_link = self._folder / fasta.name
        if not path_link.is_file():
            logger.info(f'Creating symlink for input file: {path_link}')
            path_link.symlink_to(fasta)
        return path_link
