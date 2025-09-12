from pathlib import Path

from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.bcftools.bcftoolsbase import BcftoolsBase


class BcftoolsIndex(BcftoolsBase):
    """
    Indexes bgzip compressed VCF files and BCF files.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('bcftools index', '1.17')
        self._input_key = None

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if not any(key in self._tool_inputs for key in ('BCF', 'VCF_GZ')):
            raise InvalidToolInputError("No input file found (BCF / VCF_GZ supported)")
        if len(self._tool_inputs) != 1:
            raise InvalidToolInputError("Only one type of input is supported (VCF_GZ or BCF)")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        input_key = next(k for k in self._tool_inputs.keys())
        input_file_path = self.__symlink_input(input_key) if 'symlink_input' in self._parameters else (
            self._tool_inputs[input_key][0].path)
        self.__build_command(input_file_path)
        self._execute_command()
        self.__set_output(input_key, input_file_path)

    def __symlink_input(self, key: str) -> Path:
        """
        Create a symlink for the input. This avoids cluttering the directory of the input file. This can also avoid
        errors when there are no writing permissions on the directory of the input file.
        :param key: Input key
        :return: Path to symlink input
        """
        path_link = self._folder / self._tool_inputs[key][0].path.name
        if not path_link.is_file():
            logger.info(f'Creating symlink for input file: {path_link}')
            path_link.symlink_to(self._tool_inputs[key][0].path)
        return path_link

    def __build_command(self, input_file_path: Path) -> None:
        """
        Builds the command for this tool.
        :param input_file_path: Path to the input file
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join(self._build_options(excluded_parameters=['symlink_input'])),
            str(input_file_path)])

    def __set_output(self, key: str, input_file_path: Path) -> None:
        """
        Sets the output of this tool.
        :param input_file_path: Path to the input file symlink
        :param key: Input / output key
        :return: None
        """
        self._tool_outputs[key] = [ToolIOFile(input_file_path)]
