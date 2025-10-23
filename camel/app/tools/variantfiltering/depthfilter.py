from camel.app.core.errors import InvalidParameterError
from camel.app.tools.variantfiltering.basefilter import BaseFilter


class DepthFilter(BaseFilter):
    """
    Filters variants based on absolute depth, forward depth and reverse depth.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Variant Filter: Depth', '0.1')

    @property
    def full_name(self) -> str:
        """
        Returns the full name for this filter.
        :return: Full name
        """
        return 'Depth'

    @property
    def description(self) -> str:
        """
        Returns the description for this filter.
        :return: Description
        """
        return 'Total number of mapped reads ≥<b>{}</b> at variant position, ' \
               'including ≥<b>{}</b> forward and ≥<b>{}</b> reverse read(s)'.format(
                    self._parameters['min_depth'].value if 'min_depth' in self._parameters else 0,
                    self._parameters['min_forward_depth'].value if 'min_forward_depth' in self._parameters else 0,
                    self._parameters['min_reverse_depth'].value if 'min_reverse_depth' in self._parameters else 0
                )

    def _check_parameters(self) -> None:
        """
        Checks the command line parameters.
        :return: None
        """
        if not any([param in self._parameters for param in ('min_depth', 'min_forward_depth', 'min_reverse_depth')]):
            raise InvalidParameterError('No filtering parameter found')
        super()._check_parameters()

    def _apply_filter(self) -> None:
        """
        Applies the filtering on the variants.
        :return: None
        """
        self.__build_command()
        self._execute_command()

    def __create_exclude_string(self) -> str:
        """
        Creates the exclude string.
        :return: Exclude string
        """
        parts = []
        if 'min_depth' in self._parameters:
            parts.append('DP<{}'.format(self._parameters['min_depth'].value))
        if 'min_forward_depth' in self._parameters:
            parts.append('DP4[0]+DP4[2]<{}'.format(self._parameters['min_forward_depth'].value))
        if 'min_reverse_depth' in self._parameters:
            parts.append('DP4[1]+DP4[3]<{}'.format(self._parameters['min_reverse_depth'].value))
        return ' || '.join(parts)

    def __build_command(self) -> None:
        """
        Builds the command for this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            f"--exclude '{self.__create_exclude_string()}'",
            str(self._tool_inputs['VCF_GZ'][0].path),
            '--output-type z',
            f'--output {self.output_path}'
        ] + self._get_soft_filter_options())
