from pathlib import Path

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.reports.htmlelement import HtmlElement
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.tool import Tool
from camel.app.core.utils.vcfutils import retrieve_variants


class LofreqReporter(Tool):
    """
    This tool is used to create reports for LoFreq.
    """

    SUB_FOLDER = 'output/variant_calling'

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__("LoFreq Reporter", "0.1")
        self._section = None

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'VCF' not in self._tool_inputs:
            raise InvalidToolInputError("VCF input is required")
        if 'VAL_Sample' not in self._tool_inputs:
            raise InvalidToolInputError("Sample name is required")
        if 'reference' not in self._input_informs:
            raise InvalidToolInputError("Reference information is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        tool_versions = '2.1.3.1'
        self._section = HtmlReportSection("LoFreq Variant calling", subtitle=tool_versions)
        self.__add_reference_link()
        self.__add_output_files_table()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_output_files_table(self) -> None:
        """
        Adds a table containing the output files.
        :return: None
        """
        self._section.add_header('Output VCF', 4)
        self._section.add_paragraph('Number of variants detected.')
        self.__add_vcf_table_summary(self._tool_inputs['VCF'][0].path)

    def __add_reference_link(self) -> None:
        """
        Adds a link to the reference genome that was used.;
        :return: None
        """
        info = self._input_informs['reference']
        if info.get('url') is not None:
            reference_paragraph = HtmlElement('p', 'Reference: ')
            reference_paragraph.add_html_object(HtmlElement('a', info['name'], [('href', info['url'])]))
        else:
            reference_paragraph = HtmlElement('p', f"Reference: {info['name']}")
        self._section.add_html_object(reference_paragraph)

    def __create_vcf_download_cell(self, path: Path, suffix: str) -> HtmlTableCell:
        """
        Creates a table cell that contains a download link for the VCF file.
        :param path: Path to the file
        :param suffix: Suffix to add to the filename
        :return: Table cell
        """
        filename = f"variants-lofreq-{suffix}.vcf"
        relative_path = Path(LofreqReporter.SUB_FOLDER, filename)
        self._section.add_file(path, relative_path)
        return HtmlTableCell('Download', link=str(relative_path))

    def __add_vcf_table_summary(self, path: Path) -> None:
        """
        Parses the VCF file for summary statistics.
        :return: None
        """
        all_snps = retrieve_variants(path, types=['snp'])
        all_indels = retrieve_variants(path, types=['indel'])
        vcf_cell = self.__create_vcf_download_cell(self._tool_inputs['VCF'][0].path, 'all')

        variant_table = [
            [len(all_snps), len(all_indels), vcf_cell],
        ]
        self._section.add_table(variant_table, ['SNPs', 'indels', 'VCF file'], [('class', 'data')])
