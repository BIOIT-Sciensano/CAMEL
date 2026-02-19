from distutils.util import strtobool
from pathlib import Path

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.reports.htmlelement import HtmlElement
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.tool import Tool
from camel.app.core.utils import fileutils
from camel.app.core.utils.vcfutils import retrieve_variants


class LofreqReporter(Tool):
    """
    This tool is used to create reports for LoFreq.
    """

    SUB_FOLDER = 'output/variant_calling'
    AF_TO_REPORT = [0.05, 0.1, 0.25, 0.5, 0.75, 1]

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
        if 'BAM' not in self._tool_inputs:
            raise InvalidToolInputError("BAM file required")
        if 'mapping' not in self._input_informs:
            raise InvalidToolInputError("Mapping informs are required")
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
        self.__add_mapping_info()
        self.__add_output_files_table()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_mapping_info(self) -> None:
        """
        Adds the information about the read mapping to the report.
        :return: None
        """
        self._section.add_header('Read mapping', 4)
        table_data = [[
            f"{100 * self._input_informs['map_rate']['mapping_rate']:.2f}",
            f"{self._input_informs['depth']['median_depth']:.0f}"
        ]]
        self._section.add_table(table_data, ['Mapping rate (%)', 'Median depth'], [('class', 'data')])
        if bool(strtobool(self._parameters['export_bam'].value)) is True:
            relative_path = Path('variant_calling', 'alignment-{}.bam'.format(
                fileutils.make_valid(self._tool_inputs['VAL_Sample'][0].value)))
            self._section.add_file(self._tool_inputs['BAM'][0].path, relative_path)
            self._section.add_link_to_file('Alignment (Sorted BAM)', relative_path)
        else:
            self._section.add_text("BAM file not exported, change pipeline options to include it.")

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

    def __retrieve_vars_at_specific_af(self, var_list: list, af_list: list) -> list:
        """
        From a list of vcf record, extracts variants at specific allele frequencies.
        :param var_list: variants list
        :param af_list: list of allele frequencies categories
        :return: list of categories and count
        """
        all_vars_af = [k.INFO['AF'] for k in var_list]
        dic = dict()
        for k in range(len(af_list)):
            if k == 0:
                count = len([j for j in all_vars_af if j <= af_list[k]])
                label = f"AF <= {af_list[k]:.2f}"
            else:
                count = len([j for j in all_vars_af if af_list[k - 1] < j <= af_list[k]])
                label = f"{af_list[k - 1]:.2f} < AF <= {af_list[k]:.2f}"
            dic[label] = count
        af_variant_lists = list(dic.items())
        return af_variant_lists

    def __add_vcf_table_summary(self, path: Path) -> None:
        """
        Parses the VCF file for summary statistics.
        :return: None
        """
        all_snps = retrieve_variants(path, types=['snp'])
        all_indels = retrieve_variants(path, types=['indel'])

        vcf_cell = self.__create_vcf_download_cell(self._tool_inputs['VCF'][0].path, 'all')
        variant_table_summary = [
            [len(all_snps), len(all_indels), vcf_cell],
        ]
        variant_table_afs = self.__retrieve_vars_at_specific_af(all_snps, LofreqReporter.AF_TO_REPORT)

        self._section.add_table(variant_table_summary, ['SNPs', 'indels', 'VCF file'], [('class', 'data')])
        self._section.add_paragraph('Number of SNPs detected per allele frequency categories.')
        self._section.add_table(variant_table_afs, ['AF', 'no. of SNPs'], [('class', 'data')])
