from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.errors import InvalidParameterError
from camel.app.core.tool import Tool
from camel.app.core.io.tooliofile import ToolIOFile


class Muscle(Tool):
    """
    =============================
    Muscle 3.8.31
    =============================
    MUSCLE stands for MUltiple Sequence Comparison by Log- Expectation. MUSCLE is claimed to achieve both better average accuracy and better
    speed than ClustalW2 or T-Coffee, depending on the chosen options.


    Required inputs:
    ----------------
    'FASTA'             ToolIOFile object.

    Output:
    -------
    Output can be chosen by setting the parameters
    'FASTA' (default)   ToolIOFile object.
    'LOG' (optional)    ToolIOFile object.
    'HTML'              ToolIOFile object.
    'MSF'               ToolIOFile object.
    'CLW'               ToolIOFile object.

    Parameters:
    ___________
    diags               Find diagonals (faster for similar sequences)
    maxiters <n>        Maximum number of iterations (integer, default 16)
    maxhours <h>        Maximum time to iterate in hours (default no limit)
    html                Write output in HTML format (default FASTA)
    msf                 Write output in GCG MSF format (default FASTA)
    clw                 Write output in CLUSTALW format (default FASTA)
    clwstrict           As -clw, with 'CLUSTAL W (1.81)' header
    log <logfile>       Log to file
    loga <logfile>      Append to log file
    quiet               Do not write progress messages to stderr

    """

    def __init__(self) -> None:
        """
        Initialize tool
        :return: None
        """
        super().__init__('muscle', '3.8.31')
        self._input_string = ''

    def _execute_tool(self) -> None:
        """
        Run Muscle
        :return: None
        """
        self._set_input()
        self._build_command()
        self._execute_command()
        self._set_output()

    def _set_input(self) -> None:
        """
        Set input for Muscle
        :return: None
        """
        self._input_string = f"-in {self._tool_inputs['FASTA'][0].path} "

    def _set_output(self) -> None:
        """
        Set output for Muscle
        :return: None
        """
        if 'html' in self._parameters:
            output_key = 'HTML'
        elif 'msf' in self._parameters:
            output_key = 'MSF'
        elif 'clw' in self._parameters or 'clwstrict' in self._parameters:
            output_key = 'CLW'
        else:
            output_key = 'FASTA'
        self._tool_outputs[output_key] = [ToolIOFile(self.folder / Path(self._parameters['out'].value))]
        if 'log' in self._parameters or 'loga' in self._parameters:
            try:
                self._tool_outputs['LOG'] = [ToolIOFile(self.folder / Path(self._parameters['log'].value))]
            except KeyError:
                self._tool_outputs['LOG'] = [ToolIOFile(self.folder / Path(self._parameters['loga'].value))]

    def _build_command(self) -> None:
        """
        Build command for running Muscle
        :return: None
        """
        build_options = self._build_options()
        self._command.command = " ".join([self._tool_command, self._input_string, *build_options])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _check_parameters(self) -> None:
        """
        Checks that only a single output format and log option is specified
        :return: None
        """
        super()._check_parameters()
        if len(set(self._parameters.keys()) & {'html', 'msf', 'clw', 'clwstrict'}) > 1:
            raise InvalidParameterError('Muscle tool can only have one output format specified!')
        if 'log' in self._parameters and 'loga' in self._parameters:
            raise InvalidParameterError('Muscle tool options log and loga are mutually exclusive!')
