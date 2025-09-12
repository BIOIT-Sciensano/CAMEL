from pathlib import Path

import pandas as pd

from camel.app.components.html.htmlexpandabletable import HtmlExpandableTable
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class CheckVReporter(Tool):
    """
    Creates HTML reports for the CheckV tool.
    """

    OUTPUT_FILES = [
        {'key': 'complete_genomes', 'name': 'Complete genomes'},
        {'key': 'completeness', 'name': 'Completeness'},
        {'key': 'contamination', 'name': 'Contamination'},
        {'key': 'quality_summary', 'name': 'Quality summary'}
    ]

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('CheckV reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        for output_file in CheckVReporter.OUTPUT_FILES:
            if f"TSV_{output_file['key']}" not in self._tool_inputs:
                raise InvalidToolInputError(f"TSV_{output_file['key']} input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('CheckV')

        # Add overview output table
        for output_file in CheckVReporter.OUTPUT_FILES:
            tsv_file = self._tool_inputs[f"TSV_{output_file['key']}"][0].path
            data = pd.read_table(tsv_file)
            section.add_header(output_file['name'], 3)
            section.add_html_object(HtmlExpandableTable([
                [CheckVReporter.format_cell_value(record[key]) for key in data.columns]
                for record in data.to_dict('records')
            ], [col.replace('_', ' ') for col in data.columns]))
            section.add_line_break()
            relative_path = Path('checkv') / Path(tsv_file).name
            section.add_link_to_file('Download (TSV)', relative_path)
            section.add_file(tsv_file, relative_path)
            section.add_horizontal_line()
        self._tool_outputs['HTML'] = [ToolIOValue(section)]

    @staticmethod
    def format_cell_value(value) -> str:
        """
        Formats values for the table cells.
        :param value: Value
        :return: Formatted value
        """
        if type(value) == float:
            return f'{value:.2f}'
        return value
