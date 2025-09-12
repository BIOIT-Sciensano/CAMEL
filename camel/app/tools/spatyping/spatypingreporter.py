import datetime
from typing import Any

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.spatyping.spatyping import SpaTypingHit
from camel.app.tools.tool import Tool


class SpaTypingReporter(Tool):
    """
    This tool creates a report section for the spa typing tool.
    """

    def __init__(self):
        """
        Initializes this tool.
        """
        super().__init__('Spa typing: reporter', '0.1')
        self._section = None

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'spa_typing' not in self._input_informs:
            raise InvalidToolInputError("Spa typing informs are required.")
        if 'VAL_hits' not in self._tool_inputs:
            raise InvalidToolInputError("Hits input is required.")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._section = HtmlReportSection('<i>spa</i> typing')
        self.__add_overview_table()
        self.__add_hits_table()
        self.__add_database_info()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_overview_table(self) -> None:
        """
        Adds a table with an overview of the results.
        :return: None
        """
        repeats = self._input_informs['spa_typing']['spa_type_repeats']
        table_data = [
            ['Detected <i>spa</i> type:', self._input_informs['spa_typing']['spa_type']],
            ['Repeats:', '-'.join(str(x) for x in repeats) if repeats is not None else 'NA'],
            ['Position:', self._input_informs['spa_typing'].get('genomic_coordinates', 'NA')],
            ['Strand:', self._input_informs['spa_typing'].get('strand', 'NA')]
        ]
        self._section.add_table(table_data, table_attributes=[('class', 'information')])

    def __add_hits_table(self) -> None:
        """
        Adds the table with the hit statistics.
        :return: None
        """
        self._section.add_header('Best hits', 4)
        if len(self._tool_inputs['VAL_hits']) == 0:
            self._section.add_paragraph('No hits found.')
        else:
            table_data = []
            for hit in [x.value for x in self._tool_inputs['VAL_hits']]:
                table_data.append(SpaTypingReporter.hit_to_html_row(hit))
            self._section.add_table(table_data, SpaTypingReporter.HTML_COLUMNS, [('class', 'data')])

    def __add_database_info(self) -> None:
        """
        Adds the database information to the report.
        :return: None
        """
        self._section.add_header('Database info', level=4)
        date = datetime.datetime.fromisoformat(self._input_informs['spa_typing']['db_info']['last_update_date'])
        table_data = [
            ('Last database update', date.strftime('%d-%m-%Y')),
            ('Last database change', self._input_informs['spa_typing']['db_info'].get('last_change')),
            ('Origin', self._input_informs['spa_typing']['db_info'].get('origin')),
        ]
        for i in range(0, len(table_data)):
            if table_data[i][-1] is not None:
                continue
            table_data[i] = (table_data[i][0], 'n/a')
        self._section.add_table(table_data, ['Field', 'Value'], [('class', 'data')])

    HTML_COLUMNS = ['<i>spa</i> type', 'Repeats', 'Length', '% covered', '% identity']

    @staticmethod
    def hit_to_html_row(hit: SpaTypingHit) -> list[Any]:
        """
        Converts a hit to a HTML table row.
        :param hit: Spatyping hit
        :return: HTML table row
        """
        data = [
            hit.spa_type,
            '-'.join(str(x) for x in hit.repeats),
            hit.length,
            f'{hit.percent_covered:.2f}',
            f'{hit.percent_identity:.2f}',
        ]
        if hit.is_perfect():
            return [HtmlTableCell(d, color='green') for d in data]
        else:
            return data
