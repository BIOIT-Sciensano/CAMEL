from pathlib import Path

import pandas as pd
from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltableformatter import FormatEntry, HtmlTableFormatter

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class IntegronFinderReporter(Tool):
    """
    Tool to generate HTML output for the IntegronFinder tool.
    """

    COLUMNS: list[FormatEntry] = [
        {'key': 'ID_integron', 'title': 'Integron'},
        {'key': 'ID_replicon', 'title': 'Replicon'},
        {'key': 'pos_beg', 'title': 'Start'},
        {'key': 'pos_end', 'title': 'End'},
        {'key': 'strand', 'title': 'Strand'},
        {'key': 'type_elt', 'title': 'Element'},
        {'key': 'annotation', 'title': 'Annotation'},
        {'key': 'model', 'title': 'Model'},
        {'key': 'type', 'title': 'Type'},
        {'key': 'considered_topology', 'title': 'Topology'}
    ]

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('IntegronFinder reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'integron_finder' not in self._input_informs:
            raise InvalidToolInputError('IntegronFinder informs are required')
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError('IntegronFinder output file is required (TSV)')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('IntegronFinder', subtitle=self._input_informs['integron_finder']['_name_full'])

        # Parse input data
        try:
            data_integrons = pd.read_table(self._tool_inputs['TSV'][0].path, comment='#')
            data_integrons.fillna('-')
            data_integrons['strand'] = data_integrons['strand'].apply(lambda x: '-' if x == -1 else '+')
            section.add_table(
                HtmlTableFormatter.format_table_data(data_integrons, IntegronFinderReporter.COLUMNS),
                [c['title'] for c in IntegronFinderReporter.COLUMNS],
                [('class', 'data')])
        except pd.errors.EmptyDataError:
            section.add_paragraph('No integrons detected.')

        # Add download link
        relative_path = Path('integron_finder', f"integrons_{self._parameters['name'].value}.tsv")
        section.add_link_to_file('Download (TSV)', relative_path=relative_path)
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path)
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
