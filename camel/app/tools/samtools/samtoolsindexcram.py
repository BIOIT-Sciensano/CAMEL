import shlex
from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.errors import ToolExecutionError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.samtools.samtoolsbase import SamtoolsBase


class SamtoolsIndexCram(SamtoolsBase):
    """
    Indexes sorted CRAM files.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('samtools index', version=None)

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if 'CRAM' not in self._tool_inputs:
            raise ValueError("No CRAM input file found")
        if len(self._tool_inputs['CRAM']) != 1:
            raise ValueError("Only one CRAM input file is supported")

        if 'FASTA_REF' not in self._tool_inputs:
            raise ValueError("No FASTA_REF input file found")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        input_file_path = self.__symlink_input()
        self.__build_command(input_file_path)
        self._execute_command()
        self._check_stderr(self._command)
        self._tool_outputs['CRAI'] = [ToolIOFile(Path(f"{input_file_path}.crai"))]

    def __symlink_input(self) -> Path:
        """
        Create a symlink for the input. This avoids cluttering the directory of the input file. This can also avoid
        errors when there are no writing permissions on the directory of the input file.
        :return: Path to symlink input
        """
        if 'output_filename' in self._parameters:
            basename = self._parameters['output_filename'].value
        else:
            basename = self._tool_inputs['CRAM'][0].basename
        new_path = Path(self.folder) / basename
        if (not new_path.is_symlink()) and (new_path != self._tool_inputs['CRAM'][0].path):
            new_path.symlink_to(self._tool_inputs['CRAM'][0].path)
        return new_path

    def __build_command(self, input_file_path) -> None:
        """
        Builds the command for this tool.
        seq_cache_populate.pl: Create REF_CACHE. Used when indexing a CRAM.
        The entire compound command is wrapped in bash -c so that all steps run
        inside the same (pixi) environment and the exported env vars are visible
        to samtools index.
        :param input_file_path: Path to the input file
        :return: None
        """
        fasta_ref = self._tool_inputs['FASTA_REF'][0].path
        samtools_cmd = ' '.join([
            self._tool_command,
            *self._build_options(excluded_parameters=['output_filename']),
            str(input_file_path)
        ])
        # Build a compound shell command to make sure samtools dependency is available
        inner = ' && '.join([
            f'seq_cache_populate.pl -root ./ref/cache {fasta_ref}',
            'export REF_PATH=:',
            'export REF_CACHE=./ref/cache/%2s/%2s/%s',
            samtools_cmd
        ])
        self._command.command = f'bash -c {shlex.quote(inner)}'

    def _check_stderr(self, command: Command) -> None:
        """
        Validates the stderr.
        :param command: Command to check
        :return: None
        """
        if 'unsorted positions' in command.stderr:
            raise ToolExecutionError(self.name, 'CRAM file is not sorted.')
        super()._check_stderr(command)
