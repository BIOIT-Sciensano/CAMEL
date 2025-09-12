import logging
# Disable logging for matplotlib
logging.getLogger('matplotlib').setLevel(logging.WARNING)
import matplotlib
# Disable interactive plots
matplotlib.use('agg')
from pathlib import Path

import pandas as pd
import plotnine

from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlexpandablediv import HtmlExpandableDiv
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class ReporterSegmentDownsampling(Tool):
    """
    Creates an output report for the segment downsampling.
    """

    COLUMN_DATA = {
        'depth_median': {
            'title': 'Median depth',
            'threshold_warn': 50,
            'threshold_fail': 25,
            'fmt': lambda x: f'{int(x):,}'},
        'covered_rate': {
            'title': 'Covered rate',
            'threshold_warn': 0.95,
            'threshold_fail': 0.90,
            'fmt': lambda x: f'{100 * x:.2f}'}
    }

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Reporter: Segment downsampling', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError('Stats input file is required (TSV)')
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('Reference genome input file is required (FASTA)')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        """
        section = HtmlReportSection('Per-segment downsampling', level=3)
        data_preprocessing = pd.read_table(self._tool_inputs['TSV'][0].path, na_values=['-'], keep_default_na=False)
        self.__add_section_parameters(section)
        self.__add_section_plot(section, data_preprocessing)
        self.__add_section_stats(section, data_preprocessing)
        self.__add_section_output_files(section)
        self.__add_section_info(section)
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]

    def __add_section_parameters(self, section: HtmlReportSection) -> None:
        """
        Adds the section with the tool parameters.
        :param section: Report section
        :return: None
        """
        section.add_table([
            ['Maximum depth (per segment):', ReporterSegmentDownsampling.COLUMN_DATA['depth_median']['fmt'](
                self._parameters['max_depth'].value)],
            ['Coverage threshold (warn. / fail):', '90.00% / 95.00%'],
            ['Median depth threshold (warn. / fail):', '50 / 25']
        ], table_attributes=[('class', 'information')])

    def __add_section_plot(self, section: HtmlReportSection, data_in: pd.DataFrame) -> None:
        """
        Adds the section with the plot.
        :param data_in: Input data
        :param section: Report section
        :return: None
        """
        section.add_header('Plot', 4)
        div = HtmlExpandableDiv('plot_cov', label='plot')

        # Add plot
        path_png = self.__plot_cov_metrics(data_in)
        relative_path = Path('preprocess', path_png.name)
        section.add_file(path_png, relative_path)
        div.add_html_object(HtmlElement(
            'img', attributes=[('src', str(relative_path)), ('width', '800'), ('height', '500')]))
        section.add_html_object(div)

        # Add description
        column_data = ReporterSegmentDownsampling.COLUMN_DATA
        section.add_paragraph("""
            This plot illustrates the rate of each segment covered with at least <b>{}x</b> depth and the median depth.
            The left plot includes horizontal lines representing the 90% and 95% cutoffs (warning and fail threshold,
            respectively), while the right plot displays the maximum coverage used for per segment downsampling
            (<b>{}x</b>).""".format(
            column_data['depth_median']['fmt'](self._parameters['gap_depth_cutoff'].value),
            column_data['depth_median']['fmt'](self._parameters['max_depth'].value)))

    @staticmethod
    def get_table_cell(value: float, col_data: dict) -> HtmlTableCell:
        """
        Formats the target value for the output table.
        :param value: Value
        :param col_data: Column metadata
        :return: Table cell
        """
        formatted_value = col_data['fmt'](value)
        if value < col_data['threshold_fail']:
            return HtmlTableCell(formatted_value, color='red')
        if value < col_data['threshold_warn']:
            return HtmlTableCell(formatted_value, color='yellow')
        return HtmlTableCell(formatted_value, color='green')

    def __add_section_stats(self, section: HtmlReportSection, data: pd.DataFrame) -> None:
        """
        Adds the section with the statistics.
        :param section: Section
        :param data: Statistics data
        :return: None
        """
        section.add_header('Statistics', 4)
        header = ['Segment', 'Median depth (pre)', 'Median depth (post)', 'Covered rate (pre)', 'Covered rate (post)']
        rows = []
        for row in data.to_dict('records'):
            rows.append([
                row['segment'],
                self.get_table_cell(row['depth_median_pre'], ReporterSegmentDownsampling.COLUMN_DATA['depth_median']),
                self.get_table_cell(row['depth_median_post'], ReporterSegmentDownsampling.COLUMN_DATA['depth_median']),
                self.get_table_cell(row['covered_rate_pre'], ReporterSegmentDownsampling.COLUMN_DATA['covered_rate']),
                self.get_table_cell(row['covered_rate_post'], ReporterSegmentDownsampling.COLUMN_DATA['covered_rate'])
            ])
        section.add_table(rows, header, [('class', 'data')])

    def __add_section_output_files(self, section: HtmlReportSection) -> None:
        """
        Adds the section with download links for the output files.
        :param section: Section
        :return: None
        """
        section.add_header('Output files', 4)

        # FASTA file
        relative_path_fasta = Path('preprocess', self._tool_inputs['FASTA'][0].path.name)
        section.add_file(self._tool_inputs['FASTA'][0].path, relative_path_fasta)

        # BAM File
        relative_path_bam = Path('preprocess', f"{FileSystemHelper.make_valid(self._parameters['name'].value)}.bam")
        section.add_file(self._tool_inputs['BAM'][0].path, relative_path_bam)
        path_bam_idx = next(self._tool_inputs['BAM'][0].path.parent.glob('*.bai'))
        relative_path_bai = relative_path_bam.parent / f'{relative_path_bam.name}.bai'
        section.add_file(path_bam_idx, relative_path_bai)

        # Table
        section.add_table([
            ['Reference genome', HtmlTableCell('Download (FASTA)', link=str(relative_path_fasta))],
            ['Mapping to reference genome (downsampled)', HtmlTableCell('Download (BAM)', link=str(relative_path_bam))],
            ['Mapping to reference genome (downsampled)', HtmlTableCell(
                'Download index (BAI)', link=str(relative_path_bai))]
        ], ['File', 'Download'], [('class', 'data')])

    def __add_section_info(self, section: HtmlReportSection) -> None:
        """
        Adds the section with the analysis information.
        :param section: Section
        :return: None
        """
        section.add_header('Additional information', 4)
        max_depth = self._parameters['max_depth'].value
        section.add_paragraph(
            "During the pre-processing stage, the initial dataset undergoes downsampling to achieve an approximate "
            f"coverage of <b>{int(max_depth):,}x</b> for each segment. Firstly, the reads are aligned to the "
            "reference genome, followed by downsampling of the resulting BAM file. FASTQ files are then extracted from "
            "these alignments. This stage also involves assessing the quality of the input data by extracting "
            "alignment metrics such as the median depth and the percentage of the reference genome covered. Metrics "
            "are reported before (<i>pre</i>) and after downsampling (<i>post</i>).")

    # noinspection PyTypeChecker
    def __plot_cov_metrics(self, data_in: pd.DataFrame) -> Path:
        """
        Plots the coverage metrics.
        :param data_in: Input data
        :return: Path to plot
        """
        # Melt data
        target_cols = {
            'covered_rate_pre': 'Covered rate',
            'depth_median_pre': 'Median depth'
        }
        data_in_long = pd.melt(
            data_in, id_vars='segment', value_vars=list(target_cols.keys()), var_name='metric', value_name='value')

        # Indicators
        data_indicators = pd.DataFrame([
            {'metric': 'covered_rate_pre', 'value': 0.90, 'color': 'red'},
            {'metric': 'covered_rate_pre', 'value': 0.95, 'color': 'orange'},
            {'metric': 'depth_median_pre', 'value': int(self._parameters['max_depth'].value), 'color': 'grey'}
        ])

        # Plot data
        p = plotnine.ggplot(data=data_in_long, mapping=plotnine.aes(x='segment', y='value'))
        p += plotnine.geom_col(stat='identity', color='black', fill='#3aaa35')
        p += plotnine.geom_hline(mapping=plotnine.aes(yintercept='value', color='color'), data=data_indicators, linetype='dashed')
        p += plotnine.scale_color_manual({'red': '#c95117', 'orange': '#f29d00', 'grey': '#58595b'})
        p += plotnine.facet_wrap('~ metric', scales='free', labeller=plotnine.as_labeller(target_cols))
        p += plotnine.theme(
            legend_position='none',
            panel_background=plotnine.element_rect(fill="white", colour="black"),
            panel_grid_major=plotnine.element_line(colour="#F2F2F2"),
            panel_grid_minor=plotnine.element_line(colour="#F2F2F2"),
            panel_border=plotnine.element_rect(colour='#58595b', size=1)
        )

        # Save plot to file
        path_out = self.folder / 'coverage_plot.png'
        p.save(filename=str(path_out), format='png', width=8, height=5, dpi=300)
        logging.info(f'Plot exported to: {path_out}')
        return path_out
