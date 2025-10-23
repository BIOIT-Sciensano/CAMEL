from pathlib import Path

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.errors import InvalidToolInputError, PipelineExecutionError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class NcbiHumanReadScrubberReporter(Tool):
    """
    Reporting class for the NCBI human read scrubber tool.
    """

    def __init__(self) -> None:
        """
        Initializes the tool.
        """
        super().__init__('HRRT reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'SCRUBBER' not in self._input_informs:
            raise InvalidToolInputError("The tool informs are required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('Human read removal', subtitle=self._input_informs['SCRUBBER']['_name'])
        subject = 'read pairs' if self._parameters['input_format'].value == 'fastq_pe' else 'reads' if self._parameters['input_format'].value == 'fastq_se' else 'contigs'
        count_total = self._input_informs['SCRUBBER']['statistics']['count_total']
        count_removed = self._input_informs['SCRUBBER']['statistics']['count_removed']
        if count_total == count_removed:
            raise PipelineExecutionError(
                'All reads/contigs were removed from the input file(s) during scrubbing. If this is not expected, '
                'try disabling the human read scrubbing step.'
            )
        section.add_table([
            [f'Total {subject}', f'{count_total:,}'],
            [f'Removed {subject}', f'{count_removed:,}'],
            ['Removed %', f'{100 * count_removed / count_total:.2f}'],
        ], ['Category', 'Number'], [('class', 'data')])

        if 'REMOVED' in self._tool_inputs:
            self.__add_files(section, subject)

        # Store the output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]

    def __add_files(self, section: HtmlReportSection, subject: str) -> None:
        """
        Adds the files containing the removed read pairs/reads/contigs to the report.
        :param section: Report section
        :param subject: Either read_pairs, reads or contigs
        :return: None
        """
        if self._input_informs['SCRUBBER']['statistics']['count_removed'] != 0:
            if self._parameters['input_format'].value != 'fastq_pe':
                relative_path = Path('human_read_scrubbing', self._tool_inputs['REMOVED'][0].path.name)
                section.add_file(self._tool_inputs['REMOVED'][0].path, relative_path)
                section.add_link_to_file(f'Removed {subject}', relative_path)
            else:
                for file, orientation in zip(range(2), ('forward', 'reverse')):
                    relative_path = Path('human_read_scrubbing', self._tool_inputs['REMOVED'][file].path.name)
                    section.add_file(self._tool_inputs['REMOVED'][file].path, relative_path)
                    section.add_link_to_file(f'Removed reads ({orientation})', relative_path)
