from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.config import config
from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class Sistr(Tool):
    """
    Serovar predictions from whole-genome sequence assemblies by determination of antigen gene and
    cgMLST gene alleles using BLAST.
    """
    def __init__(self) -> None:
        """
        Initialize tool.
        :return: None
        """
        super().__init__('SISTR', version=None)

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command(f'{self._tool_command} --version')
        self._execute_command(command, is_version_cmd=True)
        return command.stdout.split(' ')[-1].strip()

    def _execute_tool(self):
        """
        Execute the tool.
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        super()._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError("FASTA input is required")
        if 'DIR' not in self._tool_inputs:
            raise InvalidToolInputError("Database input is required (DIR).")

    def __set_output(self) -> None:
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs['JSON'] = [ToolIOFile(self.folder / self._parameters['output_filename'].value)]

    def __build_command(self) -> None:
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            '--output-format json',
            '--use-full-cgmlst-db',
            f"{self._param_data['tmp_dir']['option']} {config.dir_temp}",
            '--qc',
            '-vv',
            *self._build_options(),
            str(self._tool_inputs['FASTA'][0].path)
        ])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command was executed successfully.
        We don't check "if 'error' in self.stderr.lower()" here because some small warnings are wrongfully displayed by
        the tool as errors.
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
