from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Flye(Tool):
    """
    Flye is a de novo assembler for single-molecule sequencing reads, such as those produced by PacBio and ONT.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Flye', '2.9.4')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTQ' not in self._tool_inputs:
            raise InvalidToolInputError("FASTQ input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            *self._build_options(['genome_size', 'output_directory', 'threads', 'meta', 'read_error', 'min_overlap',
                                  'keep_haplotypes', 'no_alt_contigs', 'asm_coverage']),
            str(self._tool_inputs['FASTQ'][0].path),
            *self._build_options(['nano_corr', 'nano_hq', 'nano_raw'])
        ])
        self._execute_command()
        self.__set_output()

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def __set_output(self) -> None:
        """
        Sets the output of the tool.
        :return: None
        """
        dir_out = self.folder / self._parameters['output_directory'].value
        self._tool_outputs['FASTA'] = [ToolIOFile(dir_out / "assembly.fasta")]
