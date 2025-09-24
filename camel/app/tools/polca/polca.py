from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.components.files.fastautils import FastaUtils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Polca(Tool):
    """
    Polishing assembly with short reads using polca.
    """

    def __init__(self) -> None:
        """
        Initializes Polca.
        :return: None
        """
        super().__init__('POLCA', '4.1.0')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._build_command()
        self._execute_command()
        self._set_output()

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA reference is required')
        if 'FASTQ_PE' not in self._tool_inputs:
            raise InvalidToolInputError('FASTQ_PE files are required')

        if not FastaUtils.is_indexed(self._tool_inputs['FASTA'][0].path):
            raise InvalidToolInputError('FASTA reference needs to be indexed')
        super()._check_input()

    def _build_command(self) -> None:
        """
        Builds the command to run polca.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            f"-a {self._tool_inputs['FASTA'][0]}",
            f"-1 {self._tool_inputs['FASTQ_PE'][0]}",
            f"-2 {self._tool_inputs['FASTQ_PE'][1]}",
            *self._build_options()])

    def _output_is_empty(self) -> bool:
        """
        Function to check if the VCF output is empty (no variants detected).
        :return: True if output is empty
        """
        return 'POLCA has found 0 variants' in self._command.stderr

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
        """
        fasta_output = self.folder / 'output_pypolca/pypolca_corrected.fasta'
        if self._output_is_empty():
            fasta_output.symlink_to(self._tool_inputs['FASTA'][0].path)
        self._tool_outputs['FASTA'] = [ToolIOFile(fasta_output)]
