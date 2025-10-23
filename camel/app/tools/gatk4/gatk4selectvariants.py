import re

from camel.app.core.errors import InvalidToolInputError, InvalidParameterError
from camel.app.loggers import logger
from camel.app.tools.gatk4.gatk4 import GATK4


class GATK4SelectVariants(GATK4):
    """
    Class for GATK SelectVariants function
    """

    def __init__(self) -> None:
        """
        Initialize GATK4SelectVariants
        :return: None
        """
        super().__init__('gatk4 SelectVariants', '4.1.9.0')
        self._specific_parameters = ['select-type-to-exclude', 'select-type-to-include', 'selectExpressions']
        self._required_inputs = ['VCF']
        self._output_type = 'VCF'

    def _check_parameters(self) -> None:
        """
        Check the parameters set to run this tool
        :return: None
        """
        if 'selectExpressions' in self._parameters:
            select_exp = self._parameters['selectExpressions'].value
            if not select_exp:
                raise InvalidParameterError("Selection criterion of 'selectExpressions' option is required.")
            if re.search(r'\s', select_exp):
                raise InvalidParameterError(
                    f"No space in the 'selectExpressions' option JEXL expression. e.g., 'DQ>100', not 'DQ > 100'. Expression specified {select_exp}.")
        else:
            logger.warning("GATK SelectVariants running without variant selection criterion.")

    def _check_input(self) -> None:
        """
        Check input for a tool and prepare command line parameters for input
        :return: None
        """
        super()._check_input()

        if 'VCF_concordance' in self._tool_inputs and len(self._tool_inputs['VCF_concordance']) > 1:
            raise InvalidToolInputError('SelectVariant support only ONE concordance VCF file.')

        if 'VCF_discordance' in self._tool_inputs and len(self._tool_inputs['VCF_discordance']) > 1:
            raise InvalidToolInputError('SelectVariant support only ONE discordance VCF file.')

    def _set_input(self) -> None:
        """
        Set the input specification
        :return: None
        """
        super()._set_input()
        # To extract only certain samples
        if 'SAMPLES' in self._tool_inputs:
            for sample_name in self._tool_inputs['SAMPLES']:
                self._input_string += f'--sample-name {sample_name.value} '

        # To exclude certain samples
        if 'XL_SAMPLES' in self._tool_inputs:
            for sample_name in self._tool_inputs['XL_SAMPLES']:
                self._input_string += f'--exclude-sample-name {sample_name.value} '

        # Variant comparison (with those reported in another VCF file)
        # - concordance only
        if 'VCF_concordance' in self._tool_inputs:
            self._input_string += f"--concordance {self._tool_inputs['VCF_concordance'][0].path} "
        # - discordance only
        if 'VCF_discordance' in self._tool_inputs:
            self._input_string += f"--discordance {self._tool_inputs['VCF_discordance'][0].path} "

    def _set_specific_parameters(self) -> None:
        """
        Set specific parameters that need special handling
        :return: None
        """
        if 'select-type-to-exclude' in self._parameters:
            option = self._parameters['select-type-to-exclude'].option
            value = self._parameters['select-type-to-exclude'].value
            excluded_types = value.split(",")
            for ex_type in excluded_types:
                self._option_string += f"{option} {ex_type} "

        if 'select-type-to-include' in self._parameters:
            option = self._parameters['select-type-to-include'].option
            value = self._parameters['select-type-to-include'].value
            selected_types = value.split(",")
            for se_type in selected_types:
                self._option_string += f"{option} {se_type} "

        if 'selectExpressions' in self._parameters:
            select_opt = self._parameters['selectExpressions']
            self._option_string += f"{select_opt.option} {select_opt.value} "
