from distutils.util import strtobool
from pathlib import Path

import numpy as np
import pandas as pd
import plotnine

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
    AF_TO_REPORT_AND_COLOR = {0.05: '#b2182b', 0.1: '#ef8a62', 0.25: '#fddbc7', 0.5: '#d1e5f0', 0.75: '#67a9cf', 1: '#2166ac'}
    COVERAGE_THRESHOLDS_FOR_BREADTH = [1, 5, 10, 50]

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__("LoFreq Reporter", "0.1")
        self._section = None
        self._all_variants = None
        self._coverage_table = None

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
        tool_versions = 'lofreq call 2.1.3.1'
        self._section = HtmlReportSection("LoFreq Variant calling", subtitle=tool_versions)
        self.__add_reference_link()
        self.__add_mapping_info()
        self.__add_breadth_of_coverage()
        self.__add_output_files_table()
        self.__create_coverage_variant_plot()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_reference_link(self) -> None:
        """
        Adds a link to the reference genome that was used.
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

    def __add_breadth_of_coverage(self) -> None:
        """
        Adds a table with breadth of coverage
        :return: None
        """
        self._section.add_header('Breadth of coverage', 4)
        if 'TSV' in self._tool_inputs:
            res = {}
            coverage_list = [x.strip().split() for x in open(self._tool_inputs['TSV'][0].path).readlines()]
            self._coverage_table = pd.DataFrame(coverage_list, columns=['reference', 'position', 'depth'])
            self._coverage_table['depth'] = self._coverage_table['depth'].apply(pd.to_numeric)
            for cov in LofreqReporter.COVERAGE_THRESHOLDS_FOR_BREADTH:
                res[f'{cov}X'] = \
                    f'{len(self._coverage_table[self._coverage_table["depth"] >= cov]) / len(self._coverage_table) * 100:.2f}'
            self._section.add_paragraph('Breadth of coverage at different thresholds.')
            self._section.add_table(list(res.items()), ['Coverage threshold', '% covered'], [('class', 'data')])
        else:
            self._section.add_paragraph('Assay deactivated.')

    def __add_output_files_table(self) -> None:
        """
        Adds a table containing the output files.
        :return: None
        """
        self._section.add_header('Variant statistics', 4)
        self.__add_vcf_table_summary(self._tool_inputs['VCF'][0].path)

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

    @staticmethod
    def __retrieve_vars_at_specific_af(var_list: list, af_threshold: float, af_list: list) -> list:
        """
        From a list of vcf record, extracts variants at specific allele frequencies.
        :param var_list: Variants list (SNPs and Indels)
        :param af_threshold: Allele frequency threshold used as filter
        :param af_list: list of allele frequencies categories
        :return: List of categories and count
        """
        af_variants_list = []
        if var_list:
            all_indels_af = [k.INFO['AF'] for k in var_list if k.INFO.get('INDEL', False) is True]
            all_snps_af = [k.INFO['AF'] for k in var_list if k.INFO.get('INDEL', False) is False]
            for k in range(len(af_list)):
                if af_list[k] > af_threshold:
                    if k == 0:
                        count_snp = len([snp_af for snp_af in all_snps_af if snp_af <= af_list[k]])
                        count_indel = len([indel_af for indel_af in all_indels_af if indel_af <= af_list[k]])
                        label = f"AF <= {af_list[k]:.2f}"
                    else:
                        count_snp = len([snp_af for snp_af in all_snps_af if af_list[k - 1] < snp_af <= af_list[k]])
                        count_indel = len(
                            [indel_af for indel_af in all_indels_af if af_list[k - 1] < indel_af <= af_list[k]])
                        label = f"{af_list[k - 1]:.2f} < AF <= {af_list[k]:.2f}"
                    af_variants_list.append([label, count_snp, count_indel])
        return af_variants_list

    def __add_vcf_table_summary(self, path: Path) -> None:
        """
        Parses the VCF file for summary statistics.
        :param path: Path to the VCF file
        :return: None
        """
        self._all_variants = retrieve_variants(path, types=['snp', 'indel'])
        minimum_allele_frequency = self._parameters['min_af'].value if self._parameters['min_af'] else 0
        self._all_variants = [x for x in self._all_variants if x.INFO.get('AF', 0) >= minimum_allele_frequency]
        all_indels = [x for x in self._all_variants if x.INFO.get('INDEL', False) is True]

        vcf_cell = self.__create_vcf_download_cell(self._tool_inputs['VCF'][0].path, 'all')
        variant_table_summary = [
            [len(self._all_variants) - len(all_indels), len(all_indels), vcf_cell],
        ]
        variant_table_afs = self.__retrieve_vars_at_specific_af(self._all_variants, minimum_allele_frequency,
                                                                list(LofreqReporter.AF_TO_REPORT_AND_COLOR.keys()))

        self._section.add_paragraph('Number of variants detected (min AF = {:.2f})'.format(minimum_allele_frequency))
        self._section.add_table(variant_table_summary, ['Total # SNPs', 'Total # Indels', 'VCF file'],
                                [('class', 'data')])
        self._section.add_paragraph('Number of variants detected per allele frequency categories.')
        self._section.add_table(variant_table_afs, ['AF', '# SNPs', '# Indels'], [('class', 'data')])

    def __create_coverage_variant_plot(self) -> None:
        """
        Creates the coverage variant plot.
        :return: None
        """
        # Generate the depth and variants table
        # For each AF threshold in the AF_TO_REPORT table, associate them in the depth table
        depth_table = pd.DataFrame(self._coverage_table[['position', 'depth']])
        depth_table['position'] = depth_table['position'].apply(pd.to_numeric)
        depth_table['depth'] = depth_table['depth'].apply(pd.to_numeric)
        for af in list(LofreqReporter.AF_TO_REPORT_AND_COLOR.keys())[::-1]:
            var_of_interest = [k.POS for k in self._all_variants if k.INFO['AF'] >= af]
            depth_table[f'AF={af}'] = 0
            depth_table[f'AF={af}'][depth_table['position'].isin(var_of_interest)] = 1

        # Melt the depth table on the position for plotting later
        melted_depth_table = pd.melt(depth_table, 'position', depth_table.columns)
        melted_depth_table['value'] = melted_depth_table['value'].apply(pd.to_numeric)
        sub_melted = melted_depth_table[melted_depth_table['variable'] == 'depth']
        sub_melted['position'] = sub_melted['position'].apply(pd.to_numeric)
        sub_melted['value'] = sub_melted['value'].apply(pd.to_numeric)

        p = plotnine.ggplot()
        p += plotnine.geom_line(plotnine.aes(x='position', y='value'), data=sub_melted)

        # Generate the bins to plot the amount of LFV per binned genome position
        bins = np.arange(depth_table['position'].min(), depth_table['position'].max() + 250, 250)
        depth_table_binned_positions = depth_table.copy()
        depth_table_binned_positions['position_bin'] = pd.cut(depth_table_binned_positions['position'], bins=bins)

        # For all threshold in the AF_TO_REPORT table, add to the plot with a specific color
        for af in LofreqReporter.AF_TO_REPORT_AND_COLOR.keys():
            binned_df = depth_table_binned_positions.groupby('position_bin')[f'AF={af}'].agg(['sum', 'count'])
            binned_df['proportion'] = binned_df['sum'] / binned_df['count']
            binned_df = binned_df.reset_index()
            binned_df['position'] = binned_df['position_bin'].apply(lambda x: int(x.mid))
            binned_df['x_width'] = binned_df['position_bin'].apply(lambda x: x.right - x.left)
            binned_df['position'] = binned_df['position'].apply(pd.to_numeric)
            binned_df['position'] = binned_df['position'].astype(float)
            binned_df['x_width'] = binned_df['x_width'].astype(float)
            binned_df['sum'] = binned_df['sum'].apply(pd.to_numeric)

            p += plotnine.geom_bar(plotnine.aes(x='position', y='sum'), data=binned_df, stat='identity',
                                   fill=LofreqReporter.AF_TO_REPORT_AND_COLOR[af])
        p += plotnine.scale_y_log10()
        p.draw(show=True)
        p.save(f'{self._folder}/test.png', dpi=300)
        self.__add_visualization(Path(f'{self._folder}/test.png'))

    def __add_visualization(self, image_path: Path) -> None:
        """
        Adds the visualization of the mutations.
        :param image_path: Image path
        :return: None
        """
        div = HtmlElement('div', attributes=[('class', 'border_bottom')])
        div.add_header('Visualization', 4)
        relative_path = Path(LofreqReporter.SUB_FOLDER, 'figure.png')
        self._section.add_file(image_path, relative_path)
        img = HtmlElement('img', attributes=[
            ('src', str(relative_path)), ('alt', 'visualization'), ('height', '960'), ('width', '1180')])
        div.add_html_object(img)
        table_data = [
            [HtmlElement('th', 'Legend', [('colspan', '2')])],
            [HtmlTableCell('', attributes=[('style', 'background-color: #000000')]), 'Coverage'],
            [HtmlTableCell('', attributes=[('style', 'background-color: #b2182b')]), 'AF >= 0.05'],
            [HtmlTableCell('', attributes=[('style', 'background-color: #ef8a62')]), 'AF >= 0.1'],
            [HtmlTableCell('', attributes=[('style', 'background-color: #fddbc7')]), 'AF >= 0.25'],
            [HtmlTableCell('', attributes=[('style', 'background-color: #d1e5f0')]), 'AF >= 0.5'],
            [HtmlTableCell('', attributes=[('style', 'background-color: #67a9cf')]), 'AF >= 0.75']
        ]
        div.add_table(table_data, table_attributes=[('class', 'data')])
        self._section.add_html_object(div)
