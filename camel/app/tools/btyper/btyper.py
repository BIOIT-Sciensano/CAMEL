from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.core.tool import Tool


class BTyper(Tool):
    """
    In silico taxonomic classification of Bacillus cereus group isolates using assembled genomes.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('BTyper', '3.4.0')

    def _check_input(self) -> None:
        """
        Checks whether the required input files are specified.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('No FASTA input found')
        super()._check_input()

    def _build_command(self, fasta_input: Path) -> None:
        """
        Build the command to run this tool.
        :param fasta_input: Path to the FASTA input file.
        :return: None
        """
        self._command.command = ' '.join([self._tool_command, f'--input {fasta_input}', *self._build_options()])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks command output.
        :param command: Command to check.
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _set_output(self) -> None:
        """
        Collects the tool output.
        :return: None
        """
        output_filename = f'btyper3_final_results/{self._tool_inputs["FASTA"][0].path.stem}_final_results.txt'
        self._tool_outputs['TSV'] = [
            ToolIOFile(self.folder / Path(self._parameters['output_dir'].value) / Path(output_filename))]

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Symlink the input FASTA file
        fasta_input = self._folder / Path(str(self._tool_inputs['FASTA'][0])).name
        try:
            logger.info(f'Creating symlink for input FASTA file: {fasta_input}')
            fasta_input.symlink_to(self._tool_inputs["FASTA"][0].path)
        except FileExistsError:
            logger.info(f'Symlink for input FASTA file ({fasta_input}) already exists.')

        # Building the command
        self._build_command(fasta_input)
        self._execute_command()

        # Collect output
        self._set_output()
