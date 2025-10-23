from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class ReporterAmpliconClip(Tool):
    """
    Reporting class for the amplicon clipping.
    """

    COLUMNS = [
        {'title': 'Reads mapped', 'key': 'TOTAL READS'},
        {'title': 'Both clipped', 'key': 'BOTH CLIPPED'},
        {'title': 'Forward clipped', 'key': 'FORWARD CLIPPED'},
        {'title': 'Reverse clipped', 'key': 'REVERSE CLIPPED'},
        {'title': 'Not clipped', 'key': 'NOT CLIPPED'}
    ]

    def __init__(self) -> None:
        """
        Initializes the tool.
        """
        super().__init__('Reporter: ampliconclip', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'ampliconclip' not in self._input_informs:
            raise InvalidToolInputError('ampliconclip informs are required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Create report section
        section = HtmlReportSection('Primer removal', 3, subtitle=self._input_informs['ampliconclip']['_name'])

        # Information table
        section.add_table([
            [f"{self._input_informs['ampliconclip']['stats'][c['key']]:,}" for c in ReporterAmpliconClip.COLUMNS]],
            [c['title'] for c in ReporterAmpliconClip.COLUMNS], [('class', 'data')])

        # Set the tool outputs
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
