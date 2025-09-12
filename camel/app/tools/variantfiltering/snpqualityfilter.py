from camel.app.tools.variantfiltering.basefilter import BaseFilter


class SnpQualityFilter(BaseFilter):
    """
    Filters variants based on SNP quality.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Variant Filter: SNP Quality', '0.1')

    def _apply_filter(self) -> None:
        """
        Applies the filtering on the variants.
        :return: None
        """
        self.__build_command()
        self._execute_command()

    @property
    def full_name(self) -> str:
        """
        Returns the full name for this filter.
        :return: Full name
        """
        return 'SNP quality'

    @property
    def description(self) -> str:
        """
        Returns the description for this filter.
        :return: Description
        """
        return 'SNP quality ≥<b>{}</b> at variant position'.format(
            self._parameters['min_snp_quality'].value)

    def __build_command(self) -> None:
        """
        Builds the command for this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            "--exclude 'QUAL<{}'".format(self._parameters['min_snp_quality'].value),
            str(self._tool_inputs['VCF_GZ'][0].path),
            '--output-type z',
            f'--output {self.output_path}'
        ] + self._get_soft_filter_options())
