from pathlib import Path

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class ReporterTrimmingONT(Tool):
    """
    This class is used to create the trimming report for ONT reads.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('Trimming-ONT: reporter', '0.1')
        self._sub_folder = Path('read_trimming')
        self._report_section = None

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._report_section = HtmlReportSection('Read trimming', subtitle=self._input_informs['trimming']['_name'])
        self.__add_nanoplot_report('Pre-filtering', 'pre', 'HTML_PRE')
        self.__add_parameters_section()
        self.__add_statistics_section()
        self.__add_nanoplot_report('Post-filtering', 'post', 'HTML_POST')
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section)]

    def __add_nanoplot_report(self, section_title: str, suffix: str, key: str) -> None:
        """
        Adds the given NanoPlot report to the report.
        :param section_title: Title for the section
        :param suffix: Suffix for storing trimmed reads file
        :param key: Tool input key
        :return: None
        """
        self._report_section.add_header(section_title, 3)
        relative_path = self._sub_folder / f'ont_report-{suffix}.html'
        self._report_section.add_file(self._tool_inputs[key][0].path, relative_path)
        self._report_section.add_link_to_file('NanoPlot report', relative_path)

    def __add_parameters_section(self) -> None:
        """
        Adds parameters to the report.
        :return: None
        """
        self._report_section.add_header('Parameters', 3)
        self._report_section.add_table([
            ('Min. quality:', str(int(self._input_informs['trimming']['min_qual']))),
            ('Min. length:', f"{int(self._input_informs['trimming']['min_length']):,}")
        ], None, [('class', 'information')])

    def __add_statistics_section(self) -> None:
        """
        Adds the statistics section.
        :return: None
        """
        self._report_section.add_header('Statistics', 3)

        # Statistics
        report_structure = [
            {'title': 'Number of reads', 'key': 'number_of_reads', 'format': '{:,}', 'type': int},
            {'title': 'Number of bases', 'key': 'number_of_bases', 'format': '{:,.0f}', 'type': float},
            {'title': 'Mean read quality', 'key': 'mean_qual', 'format': '{:.2f}', 'type': float},
            {'title': 'Median read quality', 'key': 'median_qual', 'format': '{:.2f}', 'type': float},
            {'title': 'Mean read length', 'key': 'mean_read_length', 'format': '{:,.0f}', 'type': float},
            {'title': 'Median read length', 'key': 'median_read_length', 'format': '{:,.0f}', 'type': float},
            {'title': 'N50', 'key': 'n50', 'format': '{:,.0f}', 'type': float},
        ]
        report_data = [
            [row['title'], *[row['format'].format(row['type'](self._input_informs[step][row['key']])) for step in (
                'nanoplot_pre', 'nanoplot_post')]] for row in report_structure
        ]
        self._report_section.add_table(
            report_data, ['Metric', 'Before trimming', 'After trimming'], [('class', 'data')])
