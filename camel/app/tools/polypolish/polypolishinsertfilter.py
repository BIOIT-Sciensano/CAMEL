from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool


class PolypolishInsertFilter(Tool):
    """
    Insert size filtering aims to use read pairing to remove spurious alignments, i.e.,
    alignments of a read to the wrong instance of a repeat.

    INPUT:
    - SAM files containing read mapping information for forward and reverse reads separately

    OUTPUT:
    - Filtered SAM files
    """
    OUTPUT_PREFIX = 'alignment'

    def __init__(self) -> None:
        """
        Initializes PolypolishInsertFilter.
        :return: None
        """
        super().__init__('PolypolishInsertFilter', version=None)

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
        sam_input = [sam.path for sam in self._tool_inputs['SAM']]
        self._build_command(sam_input, PolypolishInsertFilter.OUTPUT_PREFIX)
        self._execute_command()
        self._set_output(PolypolishInsertFilter.OUTPUT_PREFIX)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid.
        :return: None
        """
        if 'SAM' not in self._tool_inputs:
            raise InvalidToolInputError('SAM alignment file is required')
        if len(self._tool_inputs['SAM']) != 2:
            raise InvalidToolInputError('Only two SAM alignment files are allowed')
        super()._check_input()

    def _build_command(self, sam_input: list, prefix: str) -> None:
        """
        Builds the command to run polypolish insert filter.
        :param sam_input: list of SAM files
        :param prefix: Prefix for the output SAM files
        :return: None
        """
        output_sam_1 = self.folder / f'{prefix}_filtered_1.sam'
        output_sam_2 = self.folder / f'{prefix}_filtered_2.sam'
        self._command.command = ' '.join([
            self._tool_command,
            *self._build_options(),
            f'--in1 {sam_input[0]}',
            f'--in2 {sam_input[1]}',
            f'--out1 {output_sam_1}',
            f'--out2 {output_sam_2}'])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _set_output(self, prefix: str) -> None:
        """
        Collects the tool output.
        :param prefix: Prefix of the output files
        :return: None
        """
        self._tool_outputs['SAM'] = [ToolIOFile(self.folder / f'{prefix}_filtered_1.sam'),
                                     ToolIOFile(self.folder / f'{prefix}_filtered_2.sam')]
