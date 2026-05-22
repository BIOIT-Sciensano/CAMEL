from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool
from camel.app.core.utils import fastautils, toolutils


class Polypolish(Tool):
    """
    Polishing assembly with short reads using polypolish.

    INPUT:
    - FASTA file of assembly
    - SAM file containing read mapping information

    OUTPUT:
    - FASTA file of polished assembly
    """
    OUTPUT_NAME = 'polished.fasta'

    def __init__(self) -> None:
        """
        Initializes Polypolish.
        :return: None
        """
        super().__init__('Polypolish', version=None)

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command(f'{self._tool_command.split()[0]} --version')
        self._execute_command(command, is_version_cmd=True)
        return command.stdout.split(' ')[-1].strip()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        fasta_input = self._tool_inputs['FASTA'][0].path
        sam_input = [sam.path for sam in self._tool_inputs['SAM']]
        fasta_output = Path(Polypolish.OUTPUT_NAME)
        self._build_command(fasta_input, sam_input, fasta_output)
        self._execute_command()
        self._set_output(fasta_output)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA reference is required')
        if 'SAM' not in self._tool_inputs:
            raise InvalidToolInputError('SAM alignment file is required')

        if len(self._tool_inputs['SAM']) > 2:
            raise InvalidToolInputError('Please input at most two SAM alignment files')
        if not fastautils.is_indexed(self._tool_inputs['FASTA'][0].path):
            raise InvalidToolInputError('FASTA reference needs to be indexed')
        super()._check_input()

    def _build_command(self, fasta_input: Path, sam_input: list[Path], fasta_output: Path) -> None:
        """
        Builds the command to run Polypolish.
        :param fasta_input: Assembly to polish
        :param sam_input: list of SAM alignment files
        :param fasta_output: Polished assembly
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            *self._build_options(),
            str(fasta_input),
            *[str(sam_file) for sam_file in sam_input],
            f'> {fasta_output}',
        ])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _set_output(self, fasta_output: Path) -> None:
        """
        Collects the tool output.
        :param fasta_output: Path of the output fasta file
        :return: None
        """
        self._tool_outputs['FASTA'] = [ToolIOFile(self.folder / f'{fasta_output}')]
