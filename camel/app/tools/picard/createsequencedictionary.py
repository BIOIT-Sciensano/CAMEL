from pathlib import Path
from camel.app.core.utils import fileutils
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.picard.picard import Picard


class CreateSequenceDictionary(Picard):
    """
    Class for Picard CreateSequenceDictionary function
    """

    def __init__(self):
        """
        Initialize a picard tool
                :return: None
        """
        super().__init__('Picard CreateSequenceDictionary', '2.23.3')
        self._required_inputs = ['FASTA_REF']
        self._specific_parameters = ['output_ext', 'symlink']
        self._fasta_file = None

    def _set_output(self) -> None:
        """
        Set the output specification
        :return: None
        """
        # NOTE that the index file should be generated under the same directory as the FASTA_REF file
        index_file_name = Path(self._fasta_file.stem + self._parameters['output_ext'].value)
        # NOTE if dictionary file exists, it should be removed otherwise Picard run will fail
        if index_file_name.exists():
            fileutils.silent_remove(index_file_name)
        self._output_string = f' O={index_file_name}'

        self._tool_outputs['FASTA_REF'] = [ToolIOFile(self._fasta_file)]

    def _set_input(self) -> None:
        """
        Function to set required and optional inputs in self._input_string
        :return: None
        """
        if 'symlink' in self._parameters:
            logger.info('Creating a symlink for the FASTA_REF input to generate SequenceDictionary locally.')
            self._fasta_file = self.__symlink_input()
        else:
            self._fasta_file = self._tool_inputs['FASTA_REF'][0].path

        logger.debug(f'symlink {self._fasta_file}')

        self._input_string += f'R={self._fasta_file} '

    def __symlink_input(self) -> Path:
        """
        Creates a symlink for the input.
        :return: Path to the symlink of the input
        """
        symlink_location = self.folder / self._tool_inputs['FASTA_REF'][0].basename
        if not symlink_location.exists():
            symlink_location.symlink_to(self._tool_inputs['FASTA_REF'][0].path)

        return symlink_location
