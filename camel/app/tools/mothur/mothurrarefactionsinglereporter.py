from pathlib import Path

from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlreportsection import HtmlReportSection

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class MothurRarefactionSingleReporter(Tool):
    """
    This tool is used to generate HTML report sections based on the Rarefaction Single (mothur) output.
    """

    TITLE = "Rarefaction (single)"

    def __init__(self) -> None:
        """
        Initializes the tool.
        """
        super().__init__('mothurrarefactionsinglereporter', '0.1')
        self._report_section = None

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'PNG' not in self._tool_inputs:
            raise InvalidToolInputError('Krona output (HTML) is required.')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the rarefaction single reporter tool.
        :return: None
        """
        self._report_section = HtmlReportSection(MothurRarefactionSingleReporter.TITLE)
        self.__add_krona_report()
        self._tool_outputs['HTML'] = [ToolIOValue(self._report_section)]

    def __add_krona_report(self) -> None:
        """
        Adds a download link to the Rarefaction figure.
        :return: None
        """
        relative_path = Path('rarefaction_single') / 'rarefaction.png'
        self._report_section.add_file(self._tool_inputs['PNG'][0].path, relative_path)
        self._report_section.add_link_to_file('Rarefaction curves (PNG)', relative_path)
