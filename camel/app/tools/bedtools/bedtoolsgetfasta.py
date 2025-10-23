from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.bedtools.bedtools import Bedtools


class BedtoolsGetFasta(Bedtools):
    """
    Bedtools GetFasta func class.
    """
    DEFAULT_OUTPUT_NAME = 'bedtools_getfata_extracted_sequences.fa'

    def __init__(self) -> None:
        """
        Initialize a samtools tool.
        :return: None
        """
        super().__init__('bedtools getfasta', '2.31.0')
        self._required_inputs = ['BED', 'FASTA']

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__set_output()
        self.__build_command()
        self._execute_command()

    def __build_command(self) -> None:
        """
        Builds the command.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join(self._build_options()),
            f'-bed {self._tool_inputs["BED"][0].path}',
            f'-fi {self._tool_inputs["FASTA"][0].path}',
            f'-fo {self.DEFAULT_OUTPUT_NAME}'
        ])

    def __set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['FASTA'] = [ToolIOFile(self._folder / self.DEFAULT_OUTPUT_NAME)]

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        self._check_required_inputs()

        if len(self._tool_inputs['BED']) != 1:
            raise InvalidToolInputError("Exactly one BED input file expected.")
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidToolInputError("Exactly one FASTA input file expected.")
        super()._check_input()
