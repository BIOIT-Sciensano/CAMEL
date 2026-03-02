import json
from distutils.util import strtobool
from pathlib import Path
from typing import Any

from camel.app.core.reports.htmlelement import HtmlElement
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool
from camel.app.core.utils import fileutils


class VariantCallingReporter(Tool):
    """
    This tool is used to create reports for the variant calling pipeline.
    """

    SUB_FOLDER = 'variant_calling'

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__("Variant Calling Reporter", "0.1")
        self._section = None

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'VCF' not in self._tool_inputs:
            raise InvalidToolInputError("VCF input is required")
        if 'VCF_filt' not in self._tool_inputs:
            raise InvalidToolInputError("Filtered VCF input is required")
        if 'VAL_Sample' not in self._tool_inputs:
            raise InvalidToolInputError("Sample name is required")
        if 'reference' not in self._input_informs:
            raise InvalidToolInputError("Reference information is required")
        if 'mapping' not in self._input_informs:
            raise InvalidToolInputError("Mapping informs are required")
        if 'JSON' not in self._tool_inputs:
            raise InvalidToolInputError("JSON output of variant filtering is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        tool_versions = ', '.join([
            self._input_informs['mapping']['_name_full'], self._input_informs['calling']['_name_full']])
        self._section = HtmlReportSection("Variant calling", subtitle=tool_versions)
        self.__add_reference_link()
        self.__add_mapping_info()
        filtering_informs = self.__parse_filtering_json(self._tool_inputs['JSON'][0].path)
        self.__add_filtering_table(filtering_informs)
        if 'BED' in self._tool_inputs:
            self.__add_regions_section(filtering_informs['region'])
        self.__add_output_files_table(filtering_informs)
        self.__add_warning()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_filtering_table(self, filtering_informs: dict[str, dict[str, Any]]) -> None:
        """
        Adds a table with the filtering information.
        :param filtering_informs: Filtering informs
        :return: None
        """
        self._section.add_header('Filtering', 4)
        filtering_table_data = []
        for filter_key in ['depth', 'snp_qual', 'mapping_qual', 'distance', 'zscore']:
            if filter_key in filtering_informs:
                stats = '{:,}/{:,}'.format(
                    filtering_informs[filter_key]['variants_out'], filtering_informs[filter_key]['variants_in'])
            else:
                stats = 'NA'
            filtering_table_data.append([
                filtering_informs[filter_key]['full_name'],
                stats,
                filtering_informs[filter_key]['description']
            ])
        self._informs['nb_variants_filtered'] = min(
            [filtering_informs[key]['variants_out'] for key in filtering_informs])
        self._section.add_table(filtering_table_data, ['Filter', 'Variants passed', 'Description'], [('class', 'data')])

    def __parse_filtering_json(self, json_path: Path) -> dict[str, dict[str, int]]:
        """
        Parses the JSON file with the filtering statistics.
        :param json_path: JSON path
        :return: Filtering stats
        """
        with json_path.open() as handle:
            return json.load(handle)

    def __add_output_files_table(self, filtering_informs: dict[str, dict[str, Any]]) -> None:
        """
        Adds a table containing the output files.
        :param filtering_informs: Filtering informs
        :return: None
        """
        self._section.add_header('Output files', 4)
        vcf = self.__create_vcf_download_cell(self._tool_inputs['VCF'][0].path, 'all')
        vcf_filtered = self.__create_vcf_download_cell(self._tool_inputs['VCF_filt'][0].path, 'filtered')
        table_data = [
            ['Unfiltered', '{:,}'.format(filtering_informs['depth']['variants_in']), vcf],
            ['Filtered (All positions)', '{:,}'.format(filtering_informs['zscore']['variants_out']), vcf_filtered],
        ]
        if 'VCF_filt_regions' in self._tool_inputs:
            vcf = self.__create_vcf_download_cell(self._tool_inputs['VCF_filt_regions'][0].path, 'filtered-regions')
            table_data.append([
                'Filtered (Regions excluded)', '{:,}'.format(filtering_informs['region']['variants_out']), vcf])
        self._section.add_table(table_data, ['', 'Number of variants', 'VCF file'], [('class', 'data')])

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

    def __create_vcf_download_cell(self, path: Path, suffix: str) -> HtmlTableCell:
        """
        Creates a table cell that contains a download link for the VCF file.
        :param path: Path to the file
        :param suffix: Suffix to add to the filename
        :return: Table cell
        """
        filename = f"variants-{self._tool_inputs['VAL_Sample'][0].value}-{suffix}.vcf"
        relative_path = Path(VariantCallingReporter.SUB_FOLDER, filename)
        self._section.add_file(path, relative_path)
        return HtmlTableCell('Download', link=str(relative_path))

    def __add_regions_section(self, region_informs: dict[str, Any]) -> None:
        """
        Adds the section with the regions that are included.
        :param region_informs: Filtering informs
        :return: None
        """
        self._section.add_header('Regions', 4)
        self._section.add_paragraph('Number of variants removed: <b>{}</b>'.format(
            region_informs['variants_in'] - region_informs['variants_out']))
        relative_path = Path(VariantCallingReporter.SUB_FOLDER, 'regions_variant_calling.bed')
        self._section.add_file(self._tool_inputs['BED'][0].path, relative_path)
        self._section.add_link_to_file('Excluded regions (BED)', relative_path)

    def __add_warning(self) -> None:
        """
        Adds a warning to the report when pseudo-reads were used as input.
        :return: None
        """
        if 'pseudo_reads' in self._parameters:
            self._section.add_warning_message("The variant calling is executed on simulated reads.")
