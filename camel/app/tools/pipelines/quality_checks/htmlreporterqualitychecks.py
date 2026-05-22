from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.tool import Tool


class HtmlReporterQualityChecks(Tool):
    """
    Tool to create HTML report for the quality checks pipeline.
    """

    COLOR_CODES = {
        'OK': 'green',
        'Warning': 'orange',
        'Skipped': 'yellow',
        'Failed': 'red'
    }

    def __init__(self) -> None:
        """
        Initialize this tool.
        :return: None
        """
        super().__init__('Quality checks reporter', '0.1')
        self._sub_folder = 'quality_control'
        self._section = HtmlReportSection('Quality checks')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__add_metric_overview()
        explanations = [
            (qc_check['full_name'], qc_check['explanation']) for qc_check in self._input_informs['qc_checks'] if
            (qc_check.get('explanation') is not None) and (qc_check.get('ori') != 'rev')]
        if len(explanations) > 0:
            self.__add_explanation_metrics(explanations)
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_metric_overview(self) -> None:
        """
        Adds the table with the metrics.
        :return: None
        """
        header = ['Metric', 'Status', 'Value', 'Warning threshold', 'Fail threshold']
        data = [[
            inf['full_name'] + (f" ({inf['ori']}.)" if inf.get('ori') is not None else ''),
            HtmlTableCell(inf['status'], color=HtmlReporterQualityChecks.COLOR_CODES[inf['status']]),
            inf['fmt_string_value'].format(inf['value']) if inf['value'] is not None else 'NA',
            inf['fmt_string_value'].format(inf['threshold_warn']),
            inf['fmt_string_value'].format(inf['threshold_fail'])
        ] for inf in self._input_informs['qc_checks']]
        self._section.add_header('Overview', 3)
        self._section.add_table(data, header, [('class', 'data')])

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        if 'qc_checks' not in self._input_informs:
            raise InvalidToolInputError("No QC checks informs found")
        super(HtmlReporterQualityChecks, self)._check_input()

    def __add_explanation_metrics(self, explanations: list[tuple[str, str]]) -> None:
        """
        Adds an explanation of the quality checks.
        :param explanations: List of explanations
        :return: None
        """
        self._section.add_header('Additional information', 3)
        self._section.add_labeled_list(explanations)
