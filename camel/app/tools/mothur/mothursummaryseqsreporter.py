from pathlib import Path

import pandas as pd
from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltablecell import HtmlTableCell

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class MothurSummarySeqsReporter(Tool):
    """
    This tool is used to generate HTML report sections based on the Mothur Summary Seqs output.
    """

    TITLE = "Mothur Summary Seqs"
    #URL_PUBMED = 'https://www.ncbi.nlm.nih.gov/pubmed/{id}'

    def __init__(self) -> None:
        """
        Initializes the tool.
        """
        super().__init__('MothurSummarySeqsReporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError('Mothur Summary Seq stats (TSV) is required.')
        if 'TSV_summary' not in self._tool_inputs:
            raise InvalidToolInputError('Mothur Summary Seq full summary (TSV) is required.')
        if 'summaryseqs' not in self._input_informs:
            raise InvalidToolInputError('Mothur Summary Seqs informs are required.')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the Mothur Summary Seqs reporter tool
        :return: None
        """
        section = HtmlReportSection(MothurSummarySeqsReporter.TITLE)

        # Add output table
        output_table = self.__parse_input_file()
        self.__add_output_table(
            section, list(output_table.columns), output_table.values.tolist(), self.get_param_value('pipeline_step'))

        # Tool output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]

        # Add link to TSV file
        relative_path = Path('summary_seqs') / self._tool_inputs['TSV_summary'][0].path.name
        section.add_link_to_file('Download complete results (TSV)', relative_path)
        section.add_file(self._tool_inputs['TSV_summary'][0].path, relative_path)

    def __parse_input_file(self) -> pd.DataFrame:
        """
        Parses the input file.
        :return: Dataframe with relevant columns
        """
        to_return = pd.read_table(
            self._tool_inputs['TSV'][0].path, header=0,
            names= ['', 'Start', 'End', 'NBases', 'Ambigs', 'Polymer', 'NumSeqs', 'empty'])
        to_return = to_return.drop('empty', axis=1).fillna('')
        return to_return

    def __add_output_table(
            self, section: HtmlReportSection, header: list[str],
            data: list[list[str | HtmlTableCell]], prefix: str) -> None:
        """
        Adds an output table to the HTML report.
        :param section: Report section
        :param header: Output table header
        :param data: Output table data
        :param prefix: Prefix for the table to add
        :return: None
        """
        if len(data) > 0:
            section.add_header(f'{prefix}', level=3)
            section.add_table(data, header, [('class', 'data')])
        else:
            section.add_paragraph('No results.')
