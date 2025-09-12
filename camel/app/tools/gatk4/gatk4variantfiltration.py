from camel.app.error import InvalidParameterError
from camel.app.tools.gatk4.gatk4 import GATK4


class GATK4VariantFiltration(GATK4):
    """
    Class for GATK VariantFiltration function
    """

    def __init__(self) -> None:
        """
        Initialize the GATK4VariantFiltration tool
        :return: None
        """
        super().__init__('gatk4 VariantFiltration', '4.1.9.0')
        self._function_name = 'VariantFiltration'
        self._specific_parameters = ['filter-names', 'filter-expressions', 'genotype-filter-names', 'genotype-filter-expressions']
        self._required_inputs = ['VCF']
        self._output_type = 'VCF'

    def _check_parameters(self) -> None:
        """
        Check the parameters set to run this tool
        :return: None
        """
        if 'filter-names' not in self._parameters and 'genotype-filter-names' not in self._parameters:
            raise InvalidParameterError(
                "GATK VariantFiltration requires at least one filter to work, both 'filter-names' and 'genotype-filter-names' specs are missing.")

        if 'filter-expressions' in self._parameters:
            if 'filter-names' not in self._parameters:
                raise InvalidParameterError(
                    "GATK VariantFiltration requires specifying names and expressions for each filter, one is missing.")
            if len(self._parameters['filter-expressions'].value.split(",")) != len(self._parameters['filter-names'].value.split(",")):
                raise InvalidParameterError(
                    "The number of provided filter names does not equal the number of expressions.")
        if 'genotype-filter-names' in self._parameters:
            if 'genotype-filter-expressions' not in self._parameters:
                raise InvalidParameterError(
                    "GATK VariantFiltration requires specifying names and expressions for each genotype filter, one is missing.")
            if len(self._parameters['genotype-filter-expressions'].value.split(",")) != len(self._parameters['genotype-filter-names'].value.split(",")):
                raise InvalidParameterError(
                    "The number of provided genotype filter names does not equal the number of expressions.")

    def _set_specific_parameters(self) -> None:
        """
        Set specific parameters that need special handling
        :return: None
        """
        if 'filter-names' in self._parameters:
            self.__set_filters(self._parameters['filter-expressions'].value.split(","), self._parameters['filter-names'].value.split(","))

        if 'genotype-filter-names' in self._parameters:
            self.__set_genotype_filters(self._parameters['genotype-filter-expressions'].value.split(","), self._parameters['genotype-filter-names'].value.split(","))

    def __set_filters(self, filter_list: list[str], filtername_list: list[str]) -> None:
        """
        Set filters (based on VCF INFO fields) to filter variants
        :param filter_list: list of filter expressions
        :param filtername_list: list of filter names
        :return: None
        """
        for expression, name in zip(filter_list, filtername_list):
            self._option_string += f'--filter-name {name} --filter-expression "{expression}" '

    def __set_genotype_filters(self, filter_list: list[str], filtername_list: list[str]) -> None:
        """
        Set genotype filters (based on VCF FORMAT fields) to filter variants
        :param filter_list: list of filter specs
        :param filtername_list: list of filter names
        :return: None
        """
        for expression, name in zip(filter_list, filtername_list):
            self._option_string += f'--genotype-filter-name {name} --genotype-filter-expression "{expression}" '
