from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class StrainGSTKmerize(Tool):
    """
    StrainGST (Strain Genome Search tool) is a tool to find close reference genomes for strains present in a sample.
    StrainGST kmerize is used to kmerize a FASTQ/A file for proper use by StrainGST run.
    """

    INPUT_KEYS = ('FASTQ', 'FASTA')

    def __init__(self) -> None:
        """
        Initializes strainGST kmerize.
        """
        super().__init__('StrainGST kmerize', '1.3.9')

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid.
        :return: None
        """
        if not any(x in self._tool_inputs for x in StrainGSTKmerize.INPUT_KEYS):
            raise InvalidToolInputError('{} input is required.'.format(' or '.join(StrainGSTKmerize.INPUT_KEYS)))
        super()._check_input()

    def _build_command(self) -> None:
        """
        Builds the command to run StrainGST kmerize.
        :return: None
        """
        input_key = 'FASTQ' if 'FASTQ' in self._tool_inputs else 'FASTA'
        self._command.command = ' '.join([
            self._tool_command,
            *self._build_options(),
            ' '.join([str(f.path) for f in self._tool_inputs[input_key]])
        ])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _set_output(self) -> None:
        """
        Collects the tool output.
        :return: None
        """
        self._tool_outputs['HDF5'] = [ToolIOFile(self.folder / self._parameters['output'].value)]

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._build_command()
        self._execute_command()
        self._set_output()
