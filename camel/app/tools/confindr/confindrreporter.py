import re
from datetime import date

from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltablecell import HtmlTableCell

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class ConFindrReporter(Tool):
    """
    Tool to generate reports for the ConFindr tool.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('ConFindr reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'confindr' not in self._input_informs:
            raise InvalidToolInputError('ConFindr informs are required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('ConFindr', subtitle=self._input_informs['confindr']['_name_full'])

        # Main table
        section.add_header('Output', 3)
        if self._input_informs['confindr']['Genus'] == 'Error processing sample':
            section.add_alert('Error processing sample, species might be missing from rMLST database', 'warning')
        elif re.match(r'^\d+$', str(self._input_informs['confindr']['NumContamSNVs'])) is None:
            section.add_alert('Error processing sample with ConFindr', 'warning')
        else:
            cell = HtmlTableCell('No', 'green') if \
                    self._input_informs['confindr']['NumContamSNVs'] < 20 else HtmlTableCell('Yes', 'red')
            section.add_table([
                [cell,
                 f"<i>{self._input_informs['confindr']['Genus']}</i>",
                 self._input_informs['confindr']['NumContamSNVs'],
                 f"{self._input_informs['confindr']['BasesExamined']:,}"
                 ]],
                ['Contaminated', 'Genus', 'Contaminated SNPs', 'Bases examined'],
                table_attributes=[('class', 'data')]
            )

        # Last update
        try:
            update_date = date(*[int(x) for x in self._input_informs['confindr']['DatabaseDownloadDate'].split('-')])
        except ValueError:
            update_date = None
        update_date_str = update_date.strftime('%d-%m-%Y') if update_date is not None else 'NA'
        section.add_paragraph(f"Last database update: {update_date_str}")
        section.add_paragraph(
            '<b>Note:</b> Samples with multiple detected genera are always classified as contaminated.')

        # add message for ONT input
        if 'input_type' in self._parameters and self._parameters['input_type'].value in ('ont', 'hybrid'):
            section.add_warning_message('ConFindr on ONT input data is experimental. Results should be interpreted with caution.')

        # Set output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
