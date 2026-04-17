from pathlib import Path

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class ProdigalReporter(Tool):
    """
    Creates an HTML report for the Prodigal analysis.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Prodigal reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA input is required')
        if 'prodigal' not in self._input_informs:
            raise InvalidToolInputError("prodigal informs are required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('Prodigal', subtitle=self._input_informs['prodigal']['_name_full'])

        # Stats
        section.add_table(
            [[f"{int(self._input_informs['prodigal']['cds'][k]):,}" for k in ('nb', 'avg_len', 'std')]],
            ['Nb. predicted CDS', 'Avg. length', 'Std. deviation'], [('class', 'data')])

        # Download link
        relative_path = Path('bacmet', 'prodigal_cds.fasta')
        section.add_file(self._tool_inputs['FASTA'][0].path, relative_path)
        section.add_link_to_file('Download (FASTA)', relative_path)

        # Tool output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
