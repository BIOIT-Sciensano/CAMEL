from camel.app.core.reports.htmlexpandablediv import HtmlExpandableDiv
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class SpoTypingReporter(Tool):
    """
    This tool creates a report for the SpoTyping output.
    """

    TITLE = 'Spoligotyping'

    def __init__(self):
        """
        Initializes this tool.
        """
        super().__init__('SpoTyping Reporter', '0.1')
        self._section = None
        self._section_spacers = HtmlReportSection('Spacer counts')

    def _check_input(self):
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'VAL_type_binary' not in self._tool_inputs:
            raise InvalidToolInputError("Binary spoligotype input (VAL_type_binary) is required")
        if 'VAL_type_octal' not in self._tool_inputs:
            raise InvalidToolInputError("Octal spoligotype input (VAL_type_octal) is required")
        if 'spotyping' not in self._input_informs:
            raise InvalidToolInputError("Spoligotype metadata is required")
        super()._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self._section = HtmlReportSection(SpoTypingReporter.TITLE, subtitle=self._input_informs['spotyping']['_name'])
        metadata = self._input_informs['spotyping']['metadata']
        table_data = [
            ['Binary:', self._tool_inputs['VAL_type_binary'][0].value],
            ['Octal:', self._tool_inputs['VAL_type_octal'][0].value],
            ['SIT number:', metadata['SIT']],
            ['SpolDB4 occurence:', metadata['total']],
            ['Label:', metadata['label']]
        ]
        self._section.add_table(table_data, table_attributes=[('class', 'information')])
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

        if 'LOG' in self._tool_inputs:
            self.__add_spacer_section()
        if 'spoligo_param' in self._tool_inputs:
            self.__add_parameter_section()

    def __add_spacer_section(self):
        """
        Creates a section that contains the number of times each spacer appeared.
        :return: None
        """
        header = ['Spacer', 'Error free hits', '1-error tolerant hits', 'Detected']
        with open(self._tool_inputs['LOG'][0].path) as handle:
            table_data = []
            for line in handle.readlines()[3:]:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                code = parts[3]
                detected = 'Yes' if code == '1' else 'No'
                color = 'green' if code == '1' else 'red'
                table_data.append([parts[0], parts[1], parts[2], HtmlTableCell(detected, color)])
        div = HtmlExpandableDiv('table-spacers', 'spacer counts')
        div.add_table(table_data, header, [('class', 'data')])
        self._section.add_html_object(div)

    @staticmethod
    def generate_empty_section() -> HtmlReportSection:
        """
        Returns a report used when this analysis is disabled.
        :return: Report section
        """
        section = HtmlReportSection(SpoTypingReporter.TITLE)
        section.add_paragraph('Analysis disabled')
        return section

    def __add_parameter_section(self) -> None:
        """
        Adds a section with the parameters used for SpoTyping.
        :return: None
        """
        self._section.add_header('Parameters', 4)
        table_data = [
            ['Min. number strict hits', self._input_informs['spoligo_param']['min_strict']],
            ['Min. number relaxed hits', self._input_informs['spoligo_param']['min_relaxed']],
            ['Downsample factor', self._input_informs['spoligo_param']['downsample_factor']]
        ]
        self._section.add_table(table_data, ['Name', 'Value'], [('class', 'data')])
