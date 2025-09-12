from pathlib import Path

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class KleborateReporter(Tool):
    """
    Creates an HTML report for Kleborate.
    """

    COLUMN_MAPPING = {
        'Typing': ['ST', 'Yersiniabactin', 'YbST'],
        'Serogroup': ['wzi', 'K_locus', 'K_locus_identity', 'O_locus', 'O_locus_identity'],
        'Virulence (overview)': [
            'virulence_score', 'Colibactin', 'CbST', 'Aerobactin', 'AbST', 'Salmochelin', 'SmST', 'RmpADC', 'RmST',
            'rmpA2'],
        'Resistance': ['resistance_score']
    }

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Kleborate reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError('TSV input is required')
        if 'kleborate' not in self._input_informs:
            raise InvalidToolInputError("Kleborate informs input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('Kleborate', subtitle=self._input_informs['kleborate']['_name'])

        # Overview table
        for sub_section, columns in KleborateReporter.COLUMN_MAPPING.items():
            section.add_header(sub_section, 3)
            section.add_table([[
                self._input_informs['kleborate'][col] for col in columns]],
                [c.replace('_', ' ') for c in columns], [('class', 'data')])

        # Download link
        relative_path = Path('kleborate', 'kleborate.tsv')
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path)
        section.add_link_to_file('Download (TSV)', relative_path)

        self._tool_outputs['HTML'] = [ToolIOValue(section)]
