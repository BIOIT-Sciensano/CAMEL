from typing import Any

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class ReporterDownsampling(Tool):
    """
    Creates output reports for the down sampling workflow.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Reporter: downsampling', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'stats' not in self._input_informs:
            raise InvalidToolInputError("Statistics input is required")
        super()._check_input()

    def __add_table_params(self, section: HtmlReportSection, stats: dict[str, Any]) -> None:
        """
        Adds the table with parameters.
        :param section: HTML report section
        :param stats: Dictionary with statistics
        :return: None
        """
        section.add_table([
            ['Reference genome size:', f"{stats['size_ref_genome']:,} bases" if 'size_ref_genome' in stats else 'n/a'],
            ['Maximum coverage:', f"{stats['coverage_target']:.2f}" if 'coverage_target' in stats else 'n/a'],
        ], table_attributes=[('class', 'information')])

    def __add_table_stats(self, section: HtmlReportSection, stats: dict[str, Any], reads_are_paired: bool) -> None:
        """
        Adds the table with downsampling statistics.
        :param section: HTML report section
        :param stats: Dictionary with statistics
        :param reads_are_paired: True if reads are paired, False otherwise
        :return: None
        """
        text_downsample = 'No downsampling required' if (stats.get('downsample_factor') is None) else \
            f"{stats['downsample_factor']:.2f}"
        section.add_table([
            [f"{int(stats['total_bases']):,}", f"{int(stats['mean_read_length']):,}",
             f"{stats['coverage_estimated']:.2f}" if 'coverage_estimated' in stats else 'n/a', text_downsample]],
            ['Total bases (fwd. + rev.)' if reads_are_paired else 'Total bases', 'Read length (avg.)',
             'Estimated coverage', 'Downsample factor'],
            [('class', 'data')]
        )

    def __add_table_downsampling(self, section: HtmlReportSection, stats: dict[str, Any], reads_are_paired: bool) -> None:
        """
        Adds the table with the downsampling information.
        :param section: HTML report section
        :param stats: Dictionary with statistics
        :param reads_are_paired: True if reads are paired, False otherwise
        :return: None
        """
        if reads_are_paired:
            header = ['Read pairs in', 'Read pairs out']
            table_data = [[f"{stats['nb_read_pairs_in']:,}", f"{self._input_informs['seqtk']['reads_count'] // 2:,}"]]
        else:
            header = ['Reads in', 'Reads out']
            table_data = [[f"{stats['nb_reads_in']:,}", f"{self._input_informs['seqtk']['reads_count']:,}"]]
        section.add_header_with_subtitle('Downsampling', 3, self._input_informs['seqtk']['_name'])
        section.add_table(table_data, header, [('class', 'data')])

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Create report section
        section = HtmlReportSection('Coverage check')
        reads_are_paired = 'is_paired' in self._parameters
        stats = self._input_informs['stats']['stats']

        # Add overview tables
        self.__add_table_params(section, stats)
        self.__add_table_stats(section, stats, reads_are_paired)

        # Add downsampling table if it was performed
        if ('seqtk' in self._input_informs) and (stats['downsample_factor'] is not None):
            self.__add_table_downsampling(section, stats, reads_are_paired)
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
