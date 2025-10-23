from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils, fastautils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool


class ALE(Tool):
    """
    ALE is the Assembly Likelihood Evaluation framework that systematically evaluates the accuracy of an assembly
    in a reference-independent manner using rigorous statistical methods.
    """

    ALE_OUTPUT = "ALE.ale"

    def __init__(self) -> None:
        """
        Initializes ALE.
        :return: None
        """
        super().__init__('ALE', '2022.05.03')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        fasta_input = self._tool_inputs['FASTA'][0].path
        bam_input = self._tool_inputs['SAM'][0].path
        ale_output = Path(ALE.ALE_OUTPUT)
        self._build_command(fasta_input, bam_input, ale_output)
        self._execute_command()
        self._set_output(ale_output)
        self._parse_output(ale_output)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA reference is required')
        if 'SAM' not in self._tool_inputs:
            raise InvalidToolInputError('SAM alignment file is required')

        if not fastautils.is_indexed(self._tool_inputs['FASTA'][0].path):
            raise InvalidToolInputError('FASTA reference needs to be indexed')
        super()._check_input()

    def _build_command(self, fasta_input: Path, sam_input: Path, ale_output: Path) -> None:
        """
        Builds the command to run ALE.
        :fasta_input: Path to the assembly FASTA file
        :sam_input: Path to the SAM file of short reads mapped against the assembly
        :ale_output: Path to the ALE output file
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command, *self._build_options(), str(sam_input), str(fasta_input), str(ale_output)])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks command output.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _set_output(self, ale_output: Path) -> None:
        """
        Collects the tool output.
        :ale_output: Path to the ALE output path
        :return: None
        """
        self._tool_outputs['ALE'] = [ToolIOFile(self.folder / ale_output)]

    def _parse_output(self, ale_output: Path) -> None:
        """
        Collects the ALE score and stores it in the informs.
        :ale_output: Path to the ALE output file
        :return: None
        """
        with open(self.folder / ale_output) as handle:
            score = handle.readline().split(':')[1]
            self._informs['ale_score'] = float(score.strip())
