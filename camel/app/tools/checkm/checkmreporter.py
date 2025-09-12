from pathlib import Path

import string

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class CheckMReporter(Tool):
    """
    Creates HTML reports for the CheckM tool.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('CheckM reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError('TSV input is required')
        if 'checkm' not in self._input_informs:
            raise InvalidToolInputError('CheckM informs are required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('CheckM')
        all_columns = list(self._input_informs['checkm']['results'][0].keys())

        # Add overview output table
        columns_overview = [col for col in all_columns if col[0] not in string.digits]
        section.add_header('Overview', 3)
        section.add_table(
            [[row[key] for key in columns_overview] for row in self._input_informs['checkm']['results']],
            columns_overview, [('class', 'data')]
        )

        # Add counts
        columns_counts = ['Bin Id'] + [col for col in all_columns if col[0] in string.digits]
        section.add_header('Counts', 3)
        section.add_table(
            [[row[key] for key in columns_counts] for row in self._input_informs['checkm']['results']],
            columns_counts, [('class', 'data')]
        )

        # Add link to TSV file
        relative_path = Path('checkm') / 'output_checkm.tsv'
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path)
        section.add_link_to_file('Download all (TSV)', relative_path)

        # Set the output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
