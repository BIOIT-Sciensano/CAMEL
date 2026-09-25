import re
import shlex

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidParameterError, InvalidToolInputError
from camel.app.loggers import logger
from camel.app.tools.gatk.gatk import GATK


class GATKSelectVariants(GATK):
    """
    Selects a subset of variants from a VCF file.

    Required inputs:
    - VCF: Variants
    - FASTA_REF: Reference genome

    Optional inputs:
    - SAMPLES: Names of the samples to include (ToolIOValue)
    - XL_SAMPLES: Names of the samples to exclude (ToolIOValue)
    - VCF_concordance: Only select the variants that are also present in this VCF file
    - VCF_discordance: Only select the variants that are not present in this VCF file
    - TXT_intervals: Intervals to restrict the analysis to

    Outputs:
    - VCF: Selected variants
    """

    # Parameters that accept a comma separated list of values, passed to GATK as a repeated option
    MULTI_VALUE_PARAMETERS = ['select_type_to_exclude', 'select_type_to_include']

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('gatk SelectVariants')

    def _check_parameters(self) -> None:
        """
        Checks if the parameters are valid.
        :return: None
        """
        super()._check_parameters()
        if 'select' not in self._parameters:
            logger.warning("No variant selection criterion ('select') specified")
        elif not self.get_param_value('select'):
            raise InvalidParameterError("The 'select' parameter requires a JEXL expression")
        # The GATK wrapper script does not quote its arguments, so an expression with whitespace would be split
        elif re.search(r'\s', self.get_param_value('select')):
            raise InvalidParameterError(
                f"The 'select' JEXL expression cannot contain whitespace (e.g. 'DP>100', not 'DP > 100'): "
                f"{self.get_param_value('select')}")

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['VCF', 'FASTA_REF'])
        for key in ('VCF_concordance', 'VCF_discordance'):
            if len(self._tool_inputs.get(key, [])) > 1:
                raise InvalidToolInputError(f"Only one '{key}' input file is supported")
        super()._check_input()

    def _get_input_options(self) -> list[str]:
        """
        Returns the command line options for the input files.
        :return: Input options
        """
        options = [*super()._get_input_options(), f"-V {self._tool_inputs['VCF'][0].path}"]
        options.extend(f'--sample_name {sample.value}' for sample in self._tool_inputs.get('SAMPLES', []))
        options.extend(f'--exclude_sample_name {sample.value}' for sample in self._tool_inputs.get('XL_SAMPLES', []))
        if 'VCF_concordance' in self._tool_inputs:
            options.append(f"--concordance {self._tool_inputs['VCF_concordance'][0].path}")
        if 'VCF_discordance' in self._tool_inputs:
            options.append(f"--discordance {self._tool_inputs['VCF_discordance'][0].path}")
        return options

    def _get_parameter_options(self) -> list[str]:
        """
        Returns the command line options for the parameters.
        :return: Parameter options
        """
        options = self._build_options(excluded_parameters=['select', *GATKSelectVariants.MULTI_VALUE_PARAMETERS])
        for key in GATKSelectVariants.MULTI_VALUE_PARAMETERS:
            if key in self._parameters:
                option = self._parameters[key].option
                options.extend(f'{option} {value}' for value in str(self.get_param_value(key)).split(','))
        if 'select' in self._parameters:
            options.append(f"{self._parameters['select'].option} {shlex.quote(self.get_param_value('select'))}")
        return options

    def _set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['VCF'] = [ToolIOFile(self.folder / self.get_param_value('output'))]
