import os
from importlib.resources import files
from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.errors import ToolExecutionError, InvalidParameterError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.mega.mltreeconstruction import MLTreeConstruction
from camel.app.core.tool import Tool


class ModelSelection(Tool):
    """
    Runs MEGA model selection.
    """

    DEFAULT_OUTPUT_NAME = 'model_selection'

    def __init__(self):
        """
        Initializes this tool.
        """
        super().__init__('MEGA: Model Selection', '10.0.4')

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['FASTA'])
        super()._check_input()

    def _check_parameters(self) -> None:
        """
        Checks if the parameters are valid.
        :return: None
        """
        if self._parameters['branch_swap_filter'].value not in (
                'None', 'Very Weak', 'Weak', 'Moderate', 'Strong', 'Very Strong'):
            raise InvalidParameterError("Branch swap filter parameter value is not valid.")
        if self._parameters['missing_data_treatment'].value not in (
                'Complete deletion', 'Use all sites', 'Partial deletion'):
            raise InvalidParameterError("Missing data treatment parameter value is not valid.")
        if self._parameters['missing_data_treatment'].value == 'Partial deletion':
            if 'site_coverage_cutoff' not in self._parameters:
                raise InvalidParameterError("No site coverage cutoff given for partial deletion")
            try:
                int(self._parameters['site_coverage_cutoff'].value)
            except ValueError:
                raise InvalidParameterError("Site coverage cutoff must be an integer")
        elif 'site_coverage_cutoff' in self._parameters:
            raise InvalidParameterError("Site coverage cutoff is only applicable for 'Partial deletion'")
        super()._check_parameters()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()
        self.__analyze_output_file()

    def __build_command(self) -> None:
        """
        Builds the command line call.
        :return: None
        """
        config_file = self.__generate_config_file()
        self._command.command = ' '.join([
            self._tool_command,
            f'-d {self._tool_inputs["FASTA"][0].path}',
            f'-a {config_file}',
            f'-o {ModelSelection.DEFAULT_OUTPUT_NAME}'
        ])

    def __generate_config_file(self) -> str:
        """
        Generates the config file.
        :return: Path to config file
        """
        path_template = Path(str(files('camel').joinpath('resources/tools/mega/model_sel_ml_nucleotide_template.mao')))
        with open(path_template) as handle:
            template = handle.read()

        config_file = os.path.join(self._folder, 'config.mao')
        with open(config_file, 'w') as handle:
            handle.write(template.format(
                branch_swap_filter=self._parameters['branch_swap_filter'].value,
                missing_data_treatment=self._parameters['missing_data_treatment'].value,
                site_coverage_cutoff=self._parameters['site_coverage_cutoff'].value if
                'site_coverage_cutoff' in self._parameters else 'Not Applicable',
                threads=self._parameters['threads'].value
            ))
        return config_file

    def __set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['CSV'] = [ToolIOFile(self._folder / f'{ModelSelection.DEFAULT_OUTPUT_NAME}.csv')]
        self._tool_outputs['TXT'] = [ToolIOFile(self._folder / f'{ModelSelection.DEFAULT_OUTPUT_NAME}_summary.txt')]

    def __analyze_output_file(self) -> None:
        """
        Analyzes the output file.
        :return: None
        """
        with open(self._tool_outputs['CSV'][0].path) as handle:
            first_line = handle.readlines()[1]
            # The line has the following structure: K2+G+I, ...
            # The model is the first part (K2), +G means Gamma categories per site, +I means invariant sites
            complete_model = first_line.split(',')[0]
            self._informs['model'] = first_line.split(',')[0].split('+')[0]
            self._informs['model_full'] = MLTreeConstruction.SUBSTITUTION_MODELS[self._informs['model']]

            complete_rates = '+'.join(complete_model.split('+')[1:])
            self._informs['rates_among_sites'] = 'U' if complete_rates == '' else complete_rates
            self._informs['rates_among_sites_full'] = MLTreeConstruction.RATES_AMONG_SITES[self._informs[
                'rates_among_sites']]

    def _check_command_output(self, command: Command) -> None:
        """
        Checks the command output to see if the program executed correctly.
        :param command: Command to check
        :return: None
        """
        if 'MEGA-CC has logged the following error:' in command.stdout:
            raise ToolExecutionError(self.name, f"MEGA-CC failed to execute: {command.stdout.strip()}")
