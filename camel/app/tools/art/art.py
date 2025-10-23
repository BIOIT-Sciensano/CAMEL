from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils, fileutils
from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.core.tool import Tool


class ART(Tool):
    """
    A simulation tool to generate synthetic next-generation sequencing reads.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('ART', '2.5.8')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Build the ART command and execute it
        self.__build_command()
        self._execute_command()
        # Compress the two output files
        fileutils.gzip_compress(self.__get_output_path('1'), Path(f"{self.__get_output_path('_1')}.gz"))
        fileutils.gzip_compress(self.__get_output_path('2'), Path(f"{self.__get_output_path('_2')}.gz"))
        # Remove the uncompressed FASTQ files
        self.__remove_file(self.__get_output_path('1'))
        self.__remove_file(self.__get_output_path('2'))
        # Set the output
        self.__set_output()

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        super()._check_input()
        if 'FASTA' in self._tool_inputs:
            if len(self._tool_inputs['FASTA']) != 1:
                raise InvalidToolInputError("FASTA input requires exactly 1 file.")
        else:
            raise ValueError("FASTA input is required")

    def __build_command(self) -> None:
        """
        Builds the command.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            f"-i {self._tool_inputs['FASTA'][0].path}",
            *self._build_options()]
        )

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command was executed successfully.
        :param command: Command to check
        :return: None
        """
        if 'error' in command.stderr.lower():
            raise ToolExecutionError(self.name, f"Command execution failed (stderr: {command.stderr}).")
        toolutils.check_tool_execution(self, command, exit_code=0)

    def __get_output_path(self, suffix: str) -> Path:
        """
        Returns the path for the output file with the given suffix.
        :param suffix: suffix for the filename
        :return: Path to the generated reads
        """
        basename = self._parameters['out'].value
        return self.folder / f"{fileutils.make_valid(basename)}{suffix}.fq"

    def __set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['FASTQ_PE'] = [
            ToolIOFile(Path(f"{self.__get_output_path('_1')}.gz")), ToolIOFile(Path(f"{self.__get_output_path('_2')}.gz"))]

    @staticmethod
    def __remove_file(input_file: Path) -> None:
        """
        Removes the uncompressed FASTQ files as these are not needed for further analyses.
        :param input_file: Input file
        :return: None
        """
        logger.info(f"Removing: {input_file}")
        if not input_file.exists():
            raise FileNotFoundError(f'Output file: {input_file} does not exist')
        input_file.unlink()
