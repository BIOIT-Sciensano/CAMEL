from pathlib import Path

from Bio import SeqIO

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class Hsp65Reporter(Tool):
    """
    This class reports the output for the hsp65 based species identification for Mycobacterium.

    INPUT
    - hits (informs): List of Blast / SRST2 hits.
    - FASTA: Fasta file with the HSP65 sequence mappings

    OUTPUT
    - VAL_HTML: HTML report output of the HSP65 assay
    """

    TITLE = '<i>hsp65</i>-based differentiation'

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Mycobacterium: hsp65 reporter', '0.1')
        self._sub_folder = Path('hsp65')

    def _check_input(self) -> None:
        """
        Checks if the tool input is valid.
        :return: None
        """
        if 'FASTA_DB' not in self._tool_inputs:
            raise InvalidToolInputError("FASTA input is required")
        if 'hits' not in self._input_informs:
            raise InvalidToolInputError("Hits input is required")
        if 'columns' not in self._input_informs:
            raise InvalidToolInputError("Column name input ('columns') is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._report_section = HtmlReportSection(Hsp65Reporter.TITLE, 3)
        self.__add_output_table(self._input_informs['hits'])
        self.__add_database(self._tool_inputs['FASTA_DB'][0].path)
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section)]

    def __add_database(self, fasta_path: Path) -> None:
        """
        Adds a download link for the database to the report.
        :param fasta_path: FASTA path
        :return: None
        """
        relative_path = Path('hsp65', fasta_path.name)
        self._report_section.add_file(fasta_path, relative_path)
        self._report_section.add_link_to_file('Database (FASTA)', relative_path)

    def __get_mapping(self) -> dict[str, str]:
        """
        Returns the mapping of the sequence id to the metadata.
        :return: Mapping
        """
        mapping = {}
        with open(self._tool_inputs['FASTA'][0].path) as handle:
            seqs = SeqIO.parse(handle, 'fasta')
            for seq in seqs:
                mapping[seq.id] = ' '.join(seq.description.split(' ')[1:])
        return mapping

    def __add_output_table(self, hits: list) -> None:
        """
        Adds the output table.
        :param hits: Detected hits
        :return: None
        """
        self._report_section.add_paragraph("Detected sequences:")
        table_data = [hit.value.to_html_row(self._report_section, self._sub_folder) for hit in hits]
        self._report_section.add_table(table_data, self._input_informs['columns'], [('class', 'data')])

    @staticmethod
    def generate_empty_section() -> HtmlReportSection:
        """
        Returns an empty report for when this analysis is disabled.
        :return: Report section
        """
        section = HtmlReportSection(Hsp65Reporter.TITLE)
        section.add_paragraph('Analysis disabled')
        return section
