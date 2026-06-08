from pathlib import Path

import numpy as np
import pandas as pd
import plotnine
from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlelement import HtmlElement
from camelcore.app.reports.htmlexpandablediv import HtmlExpandableDiv
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltablecell import HtmlTableCell
from camelcore.app.utils import fileutils, vcfutils
from camelcore.app.utils.vcfutils import retrieve_variants

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool
from camel.app.toolkits.export.tsvexporter import TsvExporter


class LofreqReporter(Tool):
    """
    This tool is used to create reports for LoFreq.
    """

    SUB_FOLDER = 'output/variant_calling'
    AF_THRESHOLDS = [
        {"af": 0.05, "color": "#b2182b"},
        {"af": 0.1, "color": "#ef8a62"},
        {"af": 0.25, "color": "#fddbc7"},
        {"af": 0.5, "color": "#d1e5f0"},
        {"af": 0.75, "color": "#67a9cf"},
        {"af": 1, "color": "#2166ac"},
    ]
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
        self._sample = None

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['VCF', 'BAM', 'TSV_list'])
        if 'mapping' not in self._input_informs:
            raise InvalidToolInputError("Mapping informs are required")
        if 'reference' not in self._input_informs:
            raise InvalidToolInputError("Reference information is required")
        if 'lofreq' not in self._input_informs:
            raise InvalidToolInputError("Lofreq information is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        tool_version = self._input_informs['lofreq']['_version']
        self._section = HtmlReportSection(
            "LoFreq Variant calling", subtitle=tool_version
        )
        self.__add_reference_link()
        self.__add_mapping_info()
        self.__add_breadth_of_coverage_table()
        self.__add_summary_variants_section()
        self.__add_coverage_visualization()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_reference_link(self) -> None:
        """
        Adds a link to the reference genome that was used.
        :return: None
        """
        info = self._input_informs['reference']
        if info.get('url') is not None:
            reference_paragraph = HtmlElement('p', 'Reference: ')
            reference_paragraph.add_html_object(
                HtmlElement('a', info['name'], [('href', info['url'])])
            )
        else:
            reference_paragraph = HtmlElement('p', f"Reference: {info['name']}")
        self._section.add_html_object(reference_paragraph)

    def __add_mapping_info(self) -> None:
        """
        Adds the information about the read mapping to the report.
        :return: None
        """
        self._section.add_header('Read mapping', 4)
        table_data = [
            [
                f"{100 * self._input_informs['map_rate']['mapping_rate']:.2f}",
                f"{int(self._input_informs['depth']['median_depth']):,}",
            ]
        ]
        self._section.add_table(
            table_data, ['Mapping rate (%)', 'Median depth'], [('class', 'data')]
        )
        if self.get_param_value('export_bam') is True:
            relative_path = Path(
                'variant_calling',
                f'alignment-{fileutils.make_valid(self.get_param_value("sample_name"))}.bam',
            )
            self._section.add_file(self._tool_inputs['BAM'][0].path, relative_path)
            self._section.add_link_to_file('Alignment (Sorted BAM)', relative_path)
        else:
            self._section.add_text(
                "BAM file not exported, change pipeline options to include it."
            )

    def __add_breadth_of_coverage_table(self) -> None:
        """
        Adds a table with breadth of coverage.
        :return: None
        """
        self._section.add_header('Breadth of coverage', 4)
        if 'TSV_depth' in self._tool_inputs:
            res = {}
            coverage_list = [
                x.strip().split()
                for x in open(self._tool_inputs['TSV_depth'][0].path).readlines()
            ]
            self._coverage_table = pd.DataFrame(
                coverage_list, columns=['reference', 'position', 'depth']
            )
            self._coverage_table['depth'] = self._coverage_table['depth'].apply(
                pd.to_numeric
            )
            len_genome = len(self._coverage_table)
            for cov in LofreqReporter.COVERAGE_THRESHOLDS_FOR_BREADTH:
                cov_fraction_above_thresh = len(self._coverage_table[self._coverage_table["depth"] >= cov]) / len_genome
                res[f'{cov}X'] = f'{cov_fraction_above_thresh * 100:.2f}'
            self._section.add_paragraph('Breadth of coverage at different thresholds.')
            self._section.add_table(
                list(res.items()),
                ['Coverage threshold', '% covered'],
                [('class', 'data')],
            )
        else:
            self._section.add_paragraph('Assay deactivated.')

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

    def __retrieve_vars_at_specific_af(self) -> pd.DataFrame:
        """
        From a list of vcf record, extracts variants at specific allele frequencies.
        :return: Pandas DF of allele frequency categories and variant counts for each category
        """
        variants = vcfutils.parse_all_variants(self._tool_inputs['VCF'][0].path)
        af_and_variants_df = pd.DataFrame(
            {'AF': v.INFO.get('AF', 0), 'var_type': v.var_type} for v in variants
        )
        if af_and_variants_df.empty:
            return af_and_variants_df

        # Construct the labels
        af_labels = [
            (f"{LofreqReporter.AF_THRESHOLDS[i - 1]['af']:.2f} < " if i > 0 else '')
            + f"AF <= {LofreqReporter.AF_THRESHOLDS[i]['af']:.2f}"
            for i in range(len(LofreqReporter.AF_THRESHOLDS))
        ]

        # Assign to bins
        af_and_variants_df["AF_bin"] = pd.cut(
            af_and_variants_df["AF"],
            bins=[0] + [row['af'] for row in LofreqReporter.AF_THRESHOLDS],
            labels=af_labels,
            include_lowest=True,
        )
        dataframe_with_counts = (
            af_and_variants_df.groupby('AF_bin', observed=False)
            .var_type.value_counts()
            .unstack(fill_value=0)
        )
        dataframe_with_counts.reset_index('AF_bin', inplace=True)
        dataframe_with_counts.rename(
            columns={
                'AF_bin': 'Allele Freq. category',
                'indel': '# Indels',
                'snp': '# SNPs',
            },
            inplace=True,
        )
        return dataframe_with_counts

    def __add_complete_list_variants_table(self) -> None:
        """
        Adds the complete list of variants to the report as a div HTML object.
        :return: None
        """
        complete_table = pd.read_table(self._tool_inputs['TSV_list'][0].path, sep='\t')
        if not complete_table.empty:
            complete_table['AF'] = complete_table['AF'].map('{:.2f}'.format)
            header_complete_table = list(complete_table.columns)
            columns_to_keep_for_report = [
                'Position',
                'Type',
                'Variant',
                'Effect',
                'Gene',
                'AF',
            ]
            sub_table_for_report = complete_table[columns_to_keep_for_report]
            div = HtmlExpandableDiv('varlist', 'Complete list of variants detected.')
            div.add_table(
                data=sub_table_for_report.values.tolist(),
                column_names=list(sub_table_for_report.columns),
                table_attributes=[('class', 'data')],
            )
            self._section.add_html_object(div)
            table_path = (
                self._folder
                / f'all_variants-{fileutils.make_valid(self.get_param_value("sample_name"))}.tsv'
            )
            TsvExporter.export(
                complete_table.values.tolist(), header_complete_table, table_path
            )
            relative_path = Path(LofreqReporter.SUB_FOLDER) / table_path.name
            self._section.add_file(table_path, relative_path)
            self._section.add_link_to_file('Download (TSV)', relative_path)
            self._section.add_paragraph(
                'This table contains all the variants detected by LoFreq '
                'with their associated effect and allele frequency.'
            )
            self._section.add_paragraph(
                'The Quality value is a phred-scaled p-value describing how likely a reported SNV is a false positive. '
                'LoFreq will only report SNVs with a p-value < 5% (i.e., quality of 13) after multiple testing correction.'
            )
        else:
            self._section.add_paragraph('No variants detected.')

    def __add_summary_variants_section(self) -> None:
        """
        Parses the VCF file for summary statistics.
        :return: None
        """
        self._section.add_header('Variant statistics summary', 4)

        # First: retrieve all variants from the VCF file
        self._all_variants = retrieve_variants(self._tool_inputs['VCF'][0].path)
        minimum_allele_frequency = (
            self._parameters['min_af'].value if self._parameters['min_af'] else 0
        )
        self._all_variants = [
            var
            for var in self._all_variants
            if var.INFO.get('AF', 0) >= minimum_allele_frequency
        ]
        all_indels = [var for var in self._all_variants if var.is_indel]
        all_snps = [var for var in self._all_variants if var.var_type == 'snp']

        # Subsection: Total number of variants detected
        vcf_cell = self.__create_vcf_download_cell(
            self._tool_inputs['VCF'][0].path, 'all'
        )
        variant_table_summary = [
            [
                len(all_snps),
                len(all_indels),
                len(self._all_variants) - len(all_snps) - len(all_indels),
                vcf_cell,
            ],
        ]
        self._section.add_paragraph(
            f'Number of variants detected (min AF = {minimum_allele_frequency:.2f})'
        )
        self._section.add_table(
            variant_table_summary,
            ['Total # SNPs', 'Total # Indels', 'Other variants', 'VCF file'],
            [('class', 'data')],
        )

        # Subsection: Variants detected by AF categories
        variant_table_afs = self.__retrieve_vars_at_specific_af()
        if not variant_table_afs.empty:
            self._section.add_paragraph(
                'Number of variants detected per allele frequency categories.'
            )
            self._section.add_table(
                variant_table_afs.values.tolist(),
                list(variant_table_afs.columns),
                [('class', 'data')],
            )

        # Subsection: Complete table of variants with effect and allele frequency
        self.__add_complete_list_variants_table()
        self._section.copy_files(self._folder)

    @staticmethod
    def __bin_and_cut_table(
        df: pd.DataFrame,
        column: str,
        window_size: int,
        column_for_binning: str = 'position',
    ) -> pd.DataFrame:
        """
        Bins the table given as input based on the positions in column_for_binning. Necessary for a cleaner figure.
        :param df: pandas DataFrame
        :param column: Column to bin
        :param window_size: Window size
        :param column_for_binning: Column to base the binning on
        :return: Binned df
        """
        column_binned_renamed = f'{column_for_binning}_bin'
        bins = np.arange(
            df[column_for_binning].min(),
            df[column_for_binning].max() + window_size,
            window_size,
        )
        table_binned = df.copy()
        table_binned[column_binned_renamed] = pd.cut(
            table_binned[column_for_binning], bins=bins
        )
        final_binned_df = table_binned.groupby(column_binned_renamed, observed=False)[column].agg(
            ['sum', 'count']
        )
        final_binned_df['proportion'] = (
            final_binned_df['sum'] / final_binned_df['count']
        )
        final_binned_df = final_binned_df.reset_index()
        final_binned_df['left'] = [
            x.left for x in final_binned_df[column_binned_renamed]
        ]
        final_binned_df['right'] = [
            x.right for x in final_binned_df[column_binned_renamed]
        ]
        final_binned_df[column_for_binning] = [
            x.mid for x in final_binned_df[column_binned_renamed]
        ]
        final_binned_df['x_width'] = final_binned_df['right'] - final_binned_df['left']

        # Set the columns as numerical (again important for plotting)
        columns_to_numerical = [column_for_binning, 'sum', 'proportion', 'x_width']
        final_binned_df[columns_to_numerical] = final_binned_df[
            columns_to_numerical
        ].apply(pd.to_numeric)
        final_binned_df[column_for_binning] = final_binned_df[
            column_for_binning
        ].astype(float)
        return final_binned_df

    def __add_coverage_visualization(self) -> None:
        """
        Creates the coverage variant plot.
        :return: None
        """
        if self._coverage_table is None:
            return

        # For each AF threshold in the AF_TO_REPORT table, associate them in the depth table
        depth_table = pd.DataFrame(
            self._coverage_table[['position', 'depth']], dtype='int'
        )

        for af in [x['af'] for x in LofreqReporter.AF_THRESHOLDS][::-1]:
            var_of_interest = [
                var.POS for var in self._all_variants if var.INFO.get('AF', 0) >= af
            ]
            depth_table[f'AF={af}'] = 0
            depth_table.loc[depth_table['position'].isin(var_of_interest), f'AF={af}'] = 1

        p = plotnine.ggplot()

        # Generate the bins to plot the average depth per binned genome position
        binned_df_depth = self.__bin_and_cut_table(depth_table, 'depth', 10)
        p += plotnine.geom_line(
            plotnine.aes(x='position', y='proportion'), data=binned_df_depth
        )

        # For all threshold in the AF_TO_REPORT table, add to the plot with a specific color
        for af_info in LofreqReporter.AF_THRESHOLDS:
            binned_df_af = self.__bin_and_cut_table(
                depth_table, f'AF={af_info["af"]}', 1500
            )
            p += plotnine.geom_bar(
                plotnine.aes(x='position', y='sum'),
                data=binned_df_af,
                stat='identity',
                fill=af_info['color'],
            )
        p += plotnine.scale_y_log10()
        p += plotnine.labs(x='Position', y='Value (log-scale)')
        p.save(f'{self._folder}/figure_coverage_and_variants.png', dpi=300)
        self.__add_visualization_section(
            Path(f'{self._folder}/figure_coverage_and_variants.png')
        )

    def __add_visualization_section(self, image_path: Path) -> None:
        """
        Adds the section containing the visualization of the mutations.
        :param image_path: Image path
        :return: None
        """
        div = HtmlElement('div', attributes=[('class', 'border_bottom')])
        div.add_header('Visualization', 4)
        relative_path = Path(LofreqReporter.SUB_FOLDER, 'figure.png')
        self._section.add_file(image_path, relative_path)
        img = HtmlElement(
            'img',
            attributes=[
                ('src', str(relative_path)),
                ('alt', 'visualization'),
                ('height', '960'),
                ('width', '1180'),
            ],
        )
        div.add_html_object(img)
        rows = [
            ("#000000", "Coverage"),
            ("#b2182b", "AF >= 0.05"),
            ("#ef8a62", "AF >= 0.1"),
            ("#fddbc7", "AF >= 0.25"),
            ("#d1e5f0", "AF >= 0.5"),
            ("#67a9cf", "AF >= 0.75"),
            ("#2166ac", "AF = 1"),
        ]
        table_data = [[HtmlElement('th', 'Legend', [('colspan', '2')])]]
        for row in rows:
            table_data.append(
                [
                    HtmlTableCell(
                        '', attributes=[('style', f'background-color: {row[0]}')]
                    ),
                    row[1],
                ]
            )
        div.add_table(table_data, table_attributes=[('class', 'data')])
        self._section.add_html_object(div)
        self._section.add_paragraph(
            'This graph shows the coverage across the provided reference genome (black).'
        )
        self._section.add_paragraph(
            'Additionally, the graph shows the proportion of variants per bin with associated '
            'allele frequencies.'
        )
