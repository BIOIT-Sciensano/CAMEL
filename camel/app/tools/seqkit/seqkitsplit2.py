from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.seqkit.seqkitbase import SeqkitBase


class SeqkitSplit2(SeqkitBase):
    """
    Seqkit split2 splits sequences into files by part size or number of parts.
    """

    INPUT_KEYS = ('FASTQ', 'FASTA')

    def __init__(self) -> None:
        """
        Initializes this tool.
                :return: None
        """
        super().__init__('Seqkit split2', version=None)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if not any(x in self._tool_inputs for x in SeqkitSplit2.INPUT_KEYS):
            raise InvalidToolInputError('{} input is required.'.format(' or '.join(SeqkitSplit2.INPUT_KEYS)))
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        dir_out = (self._folder / Path(self._parameters['output_dir'].value)).absolute()
        self.__build_command(dir_out)
        self._execute_command()
        self._collect_tool_output(dir_out)

    def __build_command(self, dir_out: Path) -> None:
        """
        Builds the command line call.
        :param dir_out: Output directory
        :return: None
        """
        input_key = 'FASTQ' if 'FASTQ' in self._tool_inputs else 'FASTA'
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join([str(f.path) for f in self._tool_inputs[input_key]]),
            f' -O {dir_out}',
            *self._build_options()
        ])

    def _collect_tool_output(self, dir_out: Path) -> None:
        """
        Collects the tool output.
        :param dir_out: Output directory
        :return: None
        """
        if 'FASTA' in self._tool_inputs:
            self._tool_outputs['FASTA'] = []
            for path_fasta in sorted(dir_out.glob('*.fasta')):
                self._tool_outputs['FASTA'].append(ToolIOFile(path_fasta))
        elif 'FASTQ' in self._tool_inputs:
            self._tool_outputs['FASTQ'] = []
            for path_fastq in sorted(dir_out.glob('*.fastq')):
                self._tool_outputs['FASTQ'].append(ToolIOFile(path_fastq))

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
