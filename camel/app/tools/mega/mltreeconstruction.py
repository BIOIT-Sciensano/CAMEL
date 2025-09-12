from importlib.resources import files
from pathlib import Path

from camel.app.command.command import Command
from camel.app.error import InvalidParameterError, ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class MLTreeConstruction(Tool):
    """
    Constructs maximum likelihood trees using MEGA.
    """

    DEFAULT_OUTPUT_NAME = 'tree_construction'

    SUBSTITUTION_MODELS = {
        'JC': 'Jukes-Cantor model',
        'K2': 'Kimura 2-parameter model',
        'T92': 'Tamura 3-parameter model',
        'HKY': 'Hasegawa-Kishino-Yano model',
        'TN93': 'Tamura-Nei model',
        'GTR': 'General Time Reversible model'
    }

    HEURISTIC_METHODS = {
        'NNI': 'Nearest-Neighbor-Interchange (NNI)',
        'SPR3': 'Subtree-Pruning-Regrafting - Fast (SPR level 3)',
        'SPR5': 'Subtree-Pruning-Regrafting - Extensive (SPR level 5)'
    }

    INITIAL_TREE = {
        'NJ_BioNJ': 'Make initial tree automatically (Default - NJ/BioNJ)',
        'Max_Parsimony': 'Make initial tree automatically (Maximum Parsimony)',
        'NJ': 'Make initial tree automatically (Neighbor Joining)',
        'BioNJ': 'Make initial tree automatically (BioNJ)'
    }

    RATES_AMONG_SITES = {
        'U': 'Uniform Rates',
        'G': 'Gamma Distributed (G)',
        'I': 'Has Invariant Sites (I)',
        'G+I': 'Gamma Distributed With Invariant Sites (G+I)'
    }

    def __init__(self):
        """
        Initializes this tool.
        """
        super().__init__('MEGA: ML Tree Construction', '10.0.4')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        super()._check_input()

    def _check_parameters(self) -> None:
        """
        Checks if the provided parameters are valid.
        :return: None
        """
        if self._parameters['test_of_phylogeny'].value not in ('None', 'Bootstrap method'):
            raise InvalidParameterError("Test of phylogeny must be either 'None' or 'Bootstrap method'")
        if self._parameters['test_of_phylogeny'].value == 'Bootstrap method':
            if 'bootstrap_replications' not in self._parameters:
                raise InvalidParameterError("Number of bootstrap replications has to be specified")
            try:
                int(self._parameters['bootstrap_replications'].value)
            except ValueError:
                raise InvalidParameterError("Number of bootstrap replications has to be an integer.")
        elif 'bootstrap_replications' in self._parameters:
            raise InvalidParameterError("Number of bootstraps specified when no bootstraps are performed")
        if self._parameters['heuristic_method'].value not in MLTreeConstruction.HEURISTIC_METHODS:
            raise InvalidParameterError("Invalid heuristic method '{}' ({} supported)".format(
                self._parameters['heuristic_method'].value, ', '.join(MLTreeConstruction.HEURISTIC_METHODS.keys())))
        if self._parameters['initial_tree'].value not in MLTreeConstruction.INITIAL_TREE:
            raise InvalidParameterError("Invalid initial tree ({} supported)".format(
                ', '.join(MLTreeConstruction.INITIAL_TREE.keys())))
        if self._parameters['rates_among_sites'].value not in ('G', 'G+I') and 'gamma_categories' in self._parameters:
            raise InvalidParameterError("Gamma categories are only used when 'G' or 'G+I' rate models are used.")
        super()._check_parameters()

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
        Builds the command line call.
        :return: None
        """
        config_file = self.__generate_config_file()
        self._command.command = ' '.join([
            self._tool_command,
            f'-d {self._tool_inputs["FASTA"][0].path}',
            f'-a {config_file}',
            f'-o {MLTreeConstruction.DEFAULT_OUTPUT_NAME}'
        ])

    def __get_parameter_value(self, name: str) -> str:
        """
        Returns the value of the parameter with the given name. Returns 'Not Applicable' if the parameter is not
        specified.
        :param name: Parameter name
        :return: Parameter value of 'Not Applicable'
        """
        if name in self._parameters:
            return self._parameters[name].value
        return 'Not Applicable'

    def __generate_config_file(self) -> Path:
        """
        Generates the config file.
        :return: Path to output file
        """
        path_template = Path(str(files('camel').joinpath('resources/tools/mega/infer_ML_nucleotide_template.mao')))
        with open(path_template) as handle:
            template = handle.read()

        config_file = self._folder / 'config.mao'
        with open(config_file, 'w') as handle:
            handle.write(template.format(
                test_of_phylogeny=self.__get_parameter_value('test_of_phylogeny'),
                bootstrap_replications=self.__get_parameter_value('bootstrap_replications'),
                model=MLTreeConstruction.SUBSTITUTION_MODELS[self._parameters['model'].value],
                heuristic_method=MLTreeConstruction.HEURISTIC_METHODS[self._parameters['heuristic_method'].value],
                initial_tree=MLTreeConstruction.INITIAL_TREE[self._parameters['initial_tree'].value],
                rates_among_sites=MLTreeConstruction.RATES_AMONG_SITES[self._parameters['rates_among_sites'].value],
                branch_swap_filter=self.__get_parameter_value('branch_swap_filter'),
                gamma_categories=self.__get_parameter_value('gamma_categories'),
                missing_data_treatment=self.__get_parameter_value('missing_data_treatment'),
                site_coverage_cutoff=self.__get_parameter_value('site_coverage_cutoff'),
                threads=self.__get_parameter_value('threads')
            ))
        return config_file

    def __set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['NWK'] = [ToolIOFile(self._folder / f'{MLTreeConstruction.DEFAULT_OUTPUT_NAME}.nwk')]
        self._tool_outputs['TXT'] = [ToolIOFile(self._folder / f'{MLTreeConstruction.DEFAULT_OUTPUT_NAME}_summary.txt')]

    def _check_command_output(self, command: Command) -> None:
        """
        Checks the command output to see if the tool executed correctly.
        :param command: Commad to check
        :return: None
        """
        if 'error' in command.stdout.lower():
            raise ToolExecutionError(self.name, f"Problem generating tree: {command.stderr}")
