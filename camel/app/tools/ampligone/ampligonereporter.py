from pathlib import Path

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class AmpliGoneReporter(Tool):
    """
    Reporting class for the AmpliGone tool.
    """

    def __init__(self) -> None:
        """
        Initializes the tool.
        """
        super().__init__('AmpliGone reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'BED' not in self._tool_inputs:
            raise InvalidToolInputError('BED input is required')
        if 'ampligone' not in self._input_informs:
            raise InvalidToolInputError("AmpliGone informs are required ('ampligone')")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Combine informs into a list
        informs_combined = [self._input_informs['ampligone']] if isinstance(self._input_informs['ampligone'], dict) \
            else self._input_informs['ampligone']

        # Create report section
        section = HtmlReportSection('Primer removal', subtitle=informs_combined[0]['_name'])

        # Information table
        section.add_table([
            ['Primers file:', informs_combined[0]['fasta_primers']],
            ['Error rate:', informs_combined[0]['error_rate']],
        ])

        # Stats table
        section.add_table([[
            f"{inf['nucleotides_removed']:,}" if inf['nucleotides_removed'] is not None else '-',
            f"{inf['percentage_removed']:.2f}" if inf['percentage_removed'] is not None else '-'
        ] for inf in informs_combined], ['Nucleotides removed', '% removed'], [('class', 'data')])

        # Output and info section
        relative_path = Path('preprocess', 'primer_locations.bed')
        section.add_link_to_file('Download primer locations (BED)', relative_path)
        section.add_file(self._tool_inputs['BED'][0].path, relative_path)

        # Set the tool outputs
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
