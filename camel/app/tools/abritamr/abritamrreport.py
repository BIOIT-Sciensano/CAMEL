from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class AbriTAMRReport(Tool):

    """
    AbritAMR: abriTAMR is an AMR gene detection pipeline that runs AMRFinderPlus on a single (or list ) of given
    isolates and collates the results into a table, separating genes identified into functionally relevant groups.
    This is the report part of the abriTAMR pipeline that needs one qc file and the output from the previous step
    of the pipeline (run)
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize tool.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('AbriTAMR report', '1.0.19', camel)
        self._species = None

    def _execute_tool(self) -> None:
        """
        Executes the tool
        :return: None
        """
        self._species = self._input_informs['ABRITAMR_RUN']['species']
        self.__build_command()
        self._execute_command()
        self.__set_output()
        self._informs['_tag'] = 'REPORT'

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        super(AbriTAMRReport, self)._check_input()
        if 'TXT_MDU_QC' not in self._tool_inputs:
            raise InvalidInputSpecificationError("MDU QC file is required")
        elif not any(key in self._tool_inputs for key in ('TXT_MATCHES', 'TXT_PARTIALS')):
            raise InvalidInputSpecificationError("AbriTAMR run outputs files must be provided")

    def __set_output(self) -> None:
        """
        Sets the name of the output files.
        :return: None
        """
        self._tool_outputs['REPORT_ABRITAMR'] = [ToolIOFile(self.folder /
                                                            f"{self._parameters['output_filename'].value}_.xlsx")]

    def __build_command(self) -> None:
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        sop = 'plus' if self._species == 'Salmonella' else 'general'
        self._command.command = ' '.join([
            self._tool_command,
            '--matches', str(self._tool_inputs['TXT_MATCHES'][0]),
            '--partials', str(self._tool_inputs['TXT_PARTIALS'][0]),
            '--qc', str(self._tool_inputs['TXT_MDU_QC'][0]),
            f'--sop {sop}',
            *self._build_options()])

    def _check_command_output(self) -> None:
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if 'error' in self.stderr.lower():
            raise ToolExecutionError(f"Command execution failed (stderr: {self.stderr}).")
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Command execution failed (Exit code: {self._command.returncode})")
