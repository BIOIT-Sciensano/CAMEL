from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.tool import Tool


class AbriTAMRReport(Tool):
    """
    AbritAMR: abriTAMR is an AMR gene detection pipeline that runs AMRFinderPlus on a single (or list ) of given
    isolates and collates the results into a table, separating genes identified into functionally relevant groups.
    This is the report part of the abriTAMR pipeline that needs one qc file and the output from the previous step
    of the pipeline (run)
    """

    def __init__(self) -> None:
        """
        Initializes tool.
        :return: None
        """
        super().__init__('AbriTAMR report', '1.1.0')
        self._species = None

    def _execute_tool(self) -> None:
        """
        Executes the tool
        :return: None
        """
        self._species = self._input_informs['abritamr_run']['species']
        self.__build_command()
        self._execute_command()
        self.__set_output()
        self._informs['_tag'] = 'REPORT'

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'TXT_mdu_qc' not in self._tool_inputs:
            raise InvalidToolInputError("MDU QC file is required")
        elif not any(key in self._tool_inputs for key in ('TXT_matches', 'TXT_partials')):
            raise InvalidToolInputError("AbriTAMR run outputs files must be provided")
        super()._check_input()

    def __set_output(self) -> None:
        """
        Sets the name of the output files.
        :return: None
        """
        self._tool_outputs['REPORT_abritamr'] = [
            ToolIOFile(self.folder / f"{self.get_param_value('output_filename')}_.xlsx")]

    def __build_command(self) -> None:
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        sop = 'plus' if self._species == 'Salmonella' else 'general'
        self._command.command = ' '.join([
            self._tool_command,
            '--matches', str(self._tool_inputs['TXT_matches'][0]),
            '--partials', str(self._tool_inputs['TXT_partials'][0]),
            '--qc', str(self._tool_inputs['TXT_mdu_qc'][0]),
            f'--sop {sop}',
            *self._build_options()])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command was executed successfully.
        :param command: Command object
        :return: None
        """
        if 'error' in command.stderr.lower():
            raise ToolExecutionError(self.name, f"Command execution failed (stderr: {command.stderr}).")
        if command.exit_code != 0:
            raise ToolExecutionError(self.name, f"Command execution failed (Exit code: {command.exit_code})")
