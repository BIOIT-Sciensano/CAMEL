from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.tool import Tool
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile


class FastANI(Tool):
    """
    FastANI is developed for fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI).
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('FastANI', '1.33')

    def _check_input(self) -> None:
        """
        Checks whether the required input files are specified.
        :return: None
        """
        authorized_keys = ['FASTA_Q', 'FASTA_R', 'FASTQ_R', 'FASTQ_Q',
                           'TSV_FASTA_Q', 'TSV_FASTA_R', 'TSV_FASTQ_R', 'TSV_FASTQ_Q']

        if len([key for key in self._tool_inputs if '_Q' in key]) > 1:
            raise InvalidToolInputError('Please input at most one query (file or sequence)')
        if len([key for key in self._tool_inputs if '_R' in key]) > 1:
            raise InvalidToolInputError('Please input at most one reference (file or sequence)')

        if len([key for key in self._tool_inputs if key in authorized_keys]) != 2:
            raise InvalidToolInputError(
                'Please check your input files - maximum one query and one reference is allowed')

        super()._check_input()

    def _build_command(self) -> None:
        """
        Build the command to run tool
        :return: None
        """
        fetching_query = next(key for key in self._tool_inputs if '_Q' in key)
        query_condition = 'TSV' in fetching_query
        fetching_reference = next(key for key in self._tool_inputs if '_R' in key)
        reference_condition = 'TSV' in fetching_reference
        if not (fetching_query and fetching_reference):
            raise InvalidToolInputError(
                f'Incorrect input found: Query={fetching_query}, Ref={fetching_reference}')

        input_str_query = '--{} {}'.format(
            'queryList' if query_condition else 'query', self._tool_inputs[fetching_query][0].path)
        input_str_reference = '--{} {}'.format(
            'refList' if reference_condition else 'ref', self._tool_inputs[fetching_reference][0].path)

        self._command.command = ' '.join([
            self._tool_command, input_str_reference, input_str_query, *self._build_options()])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _set_output(self) -> None:
        """
        Collects the tool output.
        """
        output_filename = self._parameters['output_file'].value
        self._tool_outputs['TSV'] = [ToolIOFile(self.folder / Path(output_filename))]

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._build_command()
        self._execute_command()
        self._set_output()
