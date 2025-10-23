from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.samtools.samtoolsbasepipeable import SamtoolsBasePipeable


class SamtoolsMPileup(SamtoolsBasePipeable):
    """
    Produces "pileup" textual format from an alignment.
    Note: VCF output for samtools mpileup is deprecated, for VCF output, use bcftools mpileup instead.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('samtools mpileup', '1.17')

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if 'BAM' not in self._tool_inputs:
            raise InvalidToolInputError("BAM input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        path_out = self.folder / self._parameters['output_filename'].value
        self.__build_command(path_out)
        self._execute_command()
        self._tool_outputs['PILEUP'] = [ToolIOFile(path_out)]

    def __build_command(self, path_out: Path, pipe_in: bool = False, pipe_out: bool = False) -> None:
        """
        Builds the command.
        :param path_out: Output path
        :return: None
        """
        # Initialize command
        command_parts = [self._tool_command]

        # Add input
        if pipe_in:
            command_parts.append("/dev/stdin")
        else:
            command_parts.append(' '.join(str(f.path) for f in self._tool_inputs['BAM']))

        # Add optional inputs
        if 'FASTA' in self._tool_inputs:
            command_parts.append(f'--fasta-ref {self._tool_inputs["FASTA"][0].path}')
        if 'TXT_RG' in self._tool_inputs:
            command_parts.append(f'--exlude-RG {self._tool_inputs["TXT_RG"][0].path}')
        if 'TXT_POS' in self._tool_inputs:
            command_parts.append(f'--positions {self._tool_inputs["TXT_POS"][0].path}')

        # Add output
        if not pipe_out:
            command_parts.append(f'--output {path_out}')

        # Add custom options
        command_parts.extend(self._build_options())

        # Construct command
        self._command.command = ' '.join(command_parts)

    def _check_command_output(self, command: Command) -> None:
        """
        Checks the command output.
        Supersedes function in Tool class because warnings printed to stderr can cause false abort.
        """
        self._check_stderr()
        if command.exit_code != 0:
            raise ToolExecutionError(self.name, f"Error executing {self.name}: {command.stderr}")

    def _before_pipe(self, dir_, pipe_in: bool, pipe_out: bool) -> None:
        """
        Prepares the command that will be piped.
        :param dir_: Running directory
        :param pipe_in: True if tool receives piped input
        :param pipe_out: True if tool generates piped output
        :return: None
        """
        self.__build_command(self.folder / self._parameters['output_filename'].value, pipe_in, pipe_out)

    def _after_pipe(self, stderr: str, is_last_in_pipe: bool) -> None:
        """
        Performs the required steps after executing the tool as part of a pipe.
        :param stderr: Stderr for this command in the pipe
        :param is_last_in_pipe: Boolean to indicate if this is the last step in the pipe
        :return: None
        """
        if is_last_in_pipe:
            self._tool_outputs['PILEUP'] = [ToolIOFile(self.folder / self._parameters['output_filename'].value)]
