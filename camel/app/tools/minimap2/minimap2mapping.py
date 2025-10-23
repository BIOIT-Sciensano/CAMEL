from pathlib import Path
from typing import Optional

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.piping.toolpipeable import ToolPipeable


class Minimap2Mapping(ToolPipeable):
    """
    A versatile pairwise aligner for genomic and spliced nucleotide sequences.
    """

    def __init__(self) -> None:
        """
        Initializes this tool
        :return: None
        """
        super().__init__('Minimap2', '2.26')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        path_out = self.folder / self._parameters['output_filename'].value
        self.__build_command(path_out)
        self._execute_command()
        self.__set_output(path_out)

    def __build_command(self, path_out: Optional[Path], pipe_out: bool = False) -> None:
        """
        Builds the command line call.
        :param path_out: Output filename
        :param pipe_out: If true, the output is redirected to stdout
        :return: None
        """
        parts = [
            self._tool_command,
            '-ax map-ont',
            str(self._tool_inputs['FASTA'][0].path),
            str(self._tool_inputs['FASTQ'][0].path),
            ' '.join(self._build_options(excluded_parameters=['output_filename'])),
        ]
        if not pipe_out:
            parts.append(f'> {path_out}')
        self._command.command = ' '.join(parts)

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def __set_output(self, path_out: Path) -> None:
        """
        Sets the output of this tool.
        :param path_out: Output filename
        :return: None
        """
        self._tool_outputs['SAM'] = [ToolIOFile(Path(path_out))]

    def _before_pipe(self, path_pipe_in: Path, pipe_in: bool, pipe_out: bool) -> None:
        """
        Prepares the command that will be piped.
        :param path_pipe_in: Path to the input pipe
        :param pipe_in: True if tool receives piped input
        :param pipe_out: True if tool generates piped output
        :return: None
        """
        self.__build_command(None, pipe_out)

    def _after_pipe(self, stderr: str, is_last_in_pipe: bool) -> None:
        """
        Performs the required steps after executing the tool as part of a pipe.
        :param stderr: Stderr for this command in the pipe
        :param is_last_in_pipe: Boolean to indicate if this is the last step in the pipe
        :return: None
        """
        pass
