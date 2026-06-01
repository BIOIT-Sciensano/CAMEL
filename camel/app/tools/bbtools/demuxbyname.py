from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.tool import Tool


class Demuxbyname(Tool):
    """
    Demultiplexes reads into separate FASTQ files based on substrings in the read name.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('demuxbyname.sh', '38.44')

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if 'FASTQ' not in self._tool_inputs:
            raise InvalidToolInputError("No FASTQ input file found")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __build_command(self) -> None:
        """
        Builds the command for this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            f'in={self._tool_inputs["FASTQ"][0].path}',
            *self._build_options(delimiter='='),
        ])

    def __set_output(self) -> None:
        """
        Sets the tool output by globbing for files matching the output pattern.
        :return: None
        """
        glob_pattern = self.get_param_value('out').replace('%U', '*')
        output_files = sorted(self.folder.glob(glob_pattern))
        if not output_files:
            raise ToolExecutionError(self.name, f'No output files found matching: {glob_pattern}')
        self._tool_outputs['FASTQ'] = [ToolIOFile(p) for p in output_files]
