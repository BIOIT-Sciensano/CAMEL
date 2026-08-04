from pathlib import Path

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.loggers import logger
from camel.app.tools.samtools.samtoolsbase import SamtoolsBase


class SamtoolsMerge(SamtoolsBase):
    """
    Merges BAM/SAM files.

    Required inputs:
    ----------------
    BAM / SAM: At least one bam or sam file. All files should be of the same type.

    Output:
    -------
    BAM / SAM: one bam or sam file. Output filetype depends on an output file extension; bam assumed if no valid
    extension is found.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('samtools merge', version=None)

    def _check_input(self) -> None:
        """
        Checks the validity of input file types (SAM/BAM) and sets the file type accordingly.
        :return: None
        """
        if 'BAM' not in self._tool_inputs and 'SAM' not in self._tool_inputs:
            raise InvalidToolInputError("No BAM or SAM file given as input.")
        elif 'BAM' in self._tool_inputs and 'SAM' in self._tool_inputs:
            raise InvalidToolInputError("BAM and SAM files given as input; tool can only accept one type.")
        elif 'BAM' in self._tool_inputs:
            self.__input_file_type = 'BAM'
        elif 'SAM' in self._tool_inputs:
            self.__input_file_type = 'SAM'

        super(SamtoolsBase, self)._check_input()

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
            ' '.join(self._build_options(excluded_parameters=['output_filename'])),
            self._parameters['output_filename'].value,
            ' '.join([str(file.path) for file in self._tool_inputs[self.__input_file_type]])
        ])

    def __set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        output_file_type = self.__determine_output_file_type()
        output_file_path = self.folder / self._parameters['output_filename'].value
        self._tool_outputs[output_file_type] = [ToolIOFile(Path(output_file_path))]

    def __determine_output_file_type(self) -> str:
        """
        Determines the output file type to use based on output_filename.
        :return: filetype string ("BAM" or "SAM")
        """
        if self._parameters['output_filename'].value.split(".")[-1].lower() == "bam":
            filetype = "BAM"
        elif self._parameters['output_filename'].value.split(".")[-1].lower() == "sam":
            filetype = "SAM"
        else:
            filetype = "BAM"
            logger.info("Output file format not BAM or SAM. Assuming BAM and continuing.")
        return filetype

    def _check_command_output(self, command: Command) -> None:
        """
        Checks the command output.
        Supersedes function in Tool class because warnings printed to stderr can cause false abort.
        """
        self._check_stderr(command)
        toolutils.check_tool_execution(self, command, exit_code=0)
