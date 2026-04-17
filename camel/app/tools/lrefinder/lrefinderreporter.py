from typing import Union

import pandas as pd

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class LREFinderReporter(Tool):
    """
    Tool to generate HTML output for the LRE-Finder tool.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('LRE-Finder reporter', '1.0')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'lrefinder' not in self._input_informs:
            raise InvalidToolInputError('LRE-Finder informs are required')
        super()._check_input()

    @staticmethod
    def __format_value(value: Union[str, float]) -> str:
        """
        Formats a value for the output table.
        :param value: Input value
        :return: Formatted value
        """
        if type(value) is float:
            return f'{value:.2f}'
        return value

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('LRE-Finder', subtitle=self._input_informs['lrefinder']['_name_full'])
        informs = self._input_informs['lrefinder']
        section.add_paragraph(f"Detected species: <i>{informs['species']}</i>")
        data_genes = pd.DataFrame(informs['genes'])
        section.add_header('Detected genes', 3)
        section.add_table([
            [LREFinderReporter.__format_value(row[col]) for col in data_genes.columns] for
            row in data_genes.to_dict('records')], list(data_genes.columns), [('class', 'data')])
        section.add_header('Detected mutations', 3)
        data_mutations = pd.DataFrame(informs['mutations'])
        section.add_table([
            [LREFinderReporter.__format_value(row[col]) for col in data_mutations.columns] for
            row in data_mutations.to_dict('records')], list(data_mutations.columns), [('class', 'data')])
        if 'pseudo_reads' in self._parameters:
            section.add_warning_message("The tool is executed on simulated reads.")
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
