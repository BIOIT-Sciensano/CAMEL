from collections.abc import Callable
from pathlib import Path

import pandas as pd
from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlelement import HtmlElement
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltablecell import HtmlTableCell
from camelcore.app.reports.htmltableformatter import FormatEntry, HtmlTableFormatter

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.snakemake import snakemakeutils
from camel.app.core.tool import Tool


class ReporterIterativeMapping(Tool):
    """
    Class to generate reports for the iterative mapping consensus workflow.
    """

    COLS_READ_MAPPING: list[FormatEntry] = [
        {'key': 'iter', 'title': 'Iteration'},
        {'key': 'length', 'title': 'Length', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'depth_median', 'title': 'Median depth', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'depth_iqr', 'title': 'Depth IQR', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'mapping_rate', 'title': 'Mapping rate (%)', 'fmt': lambda x: f'{100 * x:.2f}'},
        {'key': 'covered_rate', 'title': 'Covered (%)', 'fmt': lambda x: f'{100 * x:.2f}'}
    ]

    COLS_LOW_DEPTH: list[FormatEntry] = [
        {'key': 'iter', 'title': 'Iteration'},
        {'key': 'length', 'title': 'Length', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'low_cov_regions', 'title': '# low depth regions', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'low_cov_total_bp', 'title': '# low depth positions', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'low_cov_perc', 'title': 'Low depth positions (%)', 'fmt': HtmlTableFormatter.FLOAT_FMT},
    ]

    REPORT_FILES = [
        {
            'rel_path': '{s}-consensus-iter_{nb}.fasta',
            'lbl': 'Consensus sequence',
            'lbl_download': 'Download (FASTA)',
            'key': 'FASTA'
        }, {
            'rel_path': '{s}-ref-iter_{nb}.fasta',
            'lbl': 'Reference sequence',
            'lbl_download': 'Download (FASTA)',
            'key': 'FASTA_ref'
        }, {
            'rel_path': '{s}-iter_{nb}.bam',
            'lbl': 'Alignment',
            'lbl_download': 'Download (BAM)',
            'key': 'BAM'
        }, {
            'rel_path': '{s}-iter_{nb}.bam.bai',
            'lbl': 'Alignment index',
            'lbl_download': 'Download (BAI)',
            'key': None
        }, {
            'rel_path': '{s}-iter_{nb}-p1-filtered.vcf',
            'lbl': 'Filtered variants (Phase 1)',
            'lbl_download': 'Download (VCF)',
            'key': 'VCF_p1'
        }, {
            'rel_path': '{s}-iter_{nb}-p2-filtered.vcf',
            'lbl': 'Filtered variants (Phase 2)',
            'lbl_download': 'Download (VCF)',
            'key': 'VCF_p2'
        }
    ]

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Reporter: Iterative mapping', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        if 'BAM' not in self._tool_inputs:
            raise InvalidToolInputError('Mapping to consensus sequence input file is required (BAM)')
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('Consensus sequence input file is required (FASTA)')
        if 'FASTA_ref' not in self._tool_inputs:
            raise InvalidToolInputError('Reference for last iteration is required (FASTA_ref)')
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError('Stats input file is required (TSV)')
        if 'TSV_seg' not in self._tool_inputs:
            raise InvalidToolInputError('Stats input file is required (TSV_seg)')
        super()._check_input()

    @staticmethod
    def __format_value(value: str, fmt: Callable | None) -> str:
        """
        Formats the given value.
        :param value: Value
        :param fmt: Formatting function
        :return: Formatted value
        """
        if value is None:
            return '-'
        if fmt is None:
            return str(value)
        return fmt(value)

    def __add_section_iterative_mapping_stats(self, section: HtmlReportSection, data_stats: pd.DataFrame) -> None:
        """
        Adds the iterative mapping statistics section.
        :param section: Report section
        :param data_stats: Statistics data
        :return: None
        """
        # Variant calling
        keys = [
            'iter',
            'phase_1-nb_snps', 'phase_1-nb_indels', 'phase_1-nb_snps_filt', 'phase_1-nb_indels_filt',
            'phase_2-nb_snps', 'phase_2-nb_indels', 'phase_2-nb_snps_filt', 'phase_2-nb_indels_filt',
        ]
        rows = [[f'{int(r[k]):,}' for k in keys] for r in data_stats.to_dict('records')]

        div = HtmlElement('div', attributes=[('class', 'border_bottom')])
        div.add_header('Statistics (All segments)', level=3)
        div.add_header('Variant calling', level=4)
        div.add_table(
            [[HtmlElement('th'), HtmlElement('th', 'Phase 1 (bcftools)', [('colspan', 4)]),
              HtmlElement('th', 'Phase 2 (Clair3)', [('colspan', 4)])],
             [HtmlElement('th'), HtmlElement('th', 'Before filtering', [('colspan', 2)]),
              HtmlElement('th', 'After filtering', [('colspan', 2)]),
              HtmlElement('th', 'Before filtering', [('colspan', 2)]),
              HtmlElement('th', 'After filtering', [('colspan', 2)])],
             [HtmlElement('th', 'Iteration'), *([HtmlElement('th', '# SNPs'), HtmlElement('th', '# indels')] * 4)]]
            + rows, None, [('class', 'data')])

        # Read mapping
        div.add_header('Read mapping', level=4)
        div.add_table(
            HtmlTableFormatter.format_table_data(data_stats, ReporterIterativeMapping.COLS_READ_MAPPING),
            [col['title'] for col in ReporterIterativeMapping.COLS_READ_MAPPING],
            [('class', 'data')]
        )

        # Per segment stats
        relative_path = Path('iterative_mapping/stats_by_segment.tsv')
        section.add_file(self._tool_inputs['TSV_seg'][0].path, relative_path)
        div.add_link_to_file('Stats per segment (TSV)', relative_path)
        section.add_html_object(div)

    def __add_section_output_files(self, section: HtmlReportSection, dir_iter_name: str) -> None:
        """
        Adds the section with download links for the output files.
        :param section: Report section
        :param dir_iter_name: Name of the directory of the last iteration
        :return: None
        """
        # Get the working directory
        nb_iter = dir_iter_name.split('_')[-1]

        # Add files
        for row in ReporterIterativeMapping.REPORT_FILES:
            if row['key'] is None:
                continue
            path_to_file = self._tool_inputs[row['key']][0].path
            section.add_file(path_to_file, Path('iterative_mapping', str(row['rel_path']).format(
                s=self._parameters['name'].value, nb=nb_iter)))

        # Manually add BAI file
        path_bam = self._tool_inputs['BAM'][0].path
        path_bai = next(path_bam.parent.glob('*.bai'))
        relative_path_bai = Path('iterative_mapping', f"{self._parameters['name'].value}-iter_{nb_iter}.bam.bai")
        section.add_file(path_bai, relative_path_bai)

        # Add table
        div = HtmlElement('div', attributes=[('class', 'border_bottom')])
        div.add_header('Output files (last iteration)', 3)
        div.add_table([
            [row['lbl'], HtmlTableCell(row['lbl_download'], link=str(
                Path('iterative_mapping', row['rel_path'].format(s=self._parameters['name'].value, nb=nb_iter))))] for
            row in ReporterIterativeMapping.REPORT_FILES
        ], ['File', 'Download'], [('class', 'data')])
        section.add_html_object(div)

    def __add_section_low_depth(self, section: HtmlReportSection, data_stats: pd.DataFrame) -> None:
        """
        Adds the section with the low-depth regions.
        :param section: Report section
        :param data_stats: Iterative mapping statistics data
        :return: None
        """
        # Create div
        div = HtmlElement('div', attributes=[('class', 'border_bottom')])
        div.add_header('Low depth regions', 3)

        # Overview table
        div.add_table([
            ['Depth cutoff:', f"{int(self._parameters['gap_depth_cutoff'].value):,}"],
            ['Positions clipped:', f"{self._input_informs['trim_fasta']['nb_clipped']:,}"],
            ['Positions masked (N):', f"{self._input_informs['trim_fasta']['nb_masked']:,}"],
        ], None, table_attributes=[('class', 'information')])

        # Add table with low coverage stats
        div.add_table(
            HtmlTableFormatter.format_table_data(data_stats, ReporterIterativeMapping.COLS_LOW_DEPTH),
            [col['title'] for col in ReporterIterativeMapping.COLS_LOW_DEPTH],
            [('class', 'data')]
        )

        # Add link to the BED file
        nb_iter = data_stats['iter'].iloc[-1]
        path_bed_io = self.folder.parents[1] / str(data_stats['dirname'].iloc[-1]) / 'phase_2-mapping' / 'bed.io'
        path_bed = snakemakeutils.load_object(path_bed_io)[0].path
        relative_path = Path('iterative_mapping', f"{self._parameters['name'].value}-iter_{nb_iter}-low_depth.bed")
        section.add_file(path_bed, relative_path)
        div.add_link_to_file('Download low coverage regions (BED)', relative_path)

        # Add explanation
        div.add_paragraph(
            """The workflow identifies regions as missing when (1) all corresponding positions have a depth of
            <b>{gap_cov}</b> or less, and (2) the length of the region is longer than <b>{gap_len}</b> bp. When missing
            regions are located at the beginning or end of a segment, they are removed or "clipped" from the consensus
            sequence. However, if the missing regions appear within the segment, they are substituted with the letter N
            to represent their absence.""".format(
                gap_cov=self._parameters['gap_depth_cutoff'].value, gap_len=self._parameters['gap_len_cutoff'].value))
        section.add_html_object(div)

    def __add_section_additional_info(self, section: HtmlReportSection) -> None:
        """
        Adds the section with additional information on the assay.
        :param section: Report section
        :return: None
        """
        div = HtmlElement('div')
        div.add_header('Additional information', 3)
        div.add_paragraph(
            """
            The iterative mapping workflow iteratively reconstructs the consensus sequence. In each iteration, the
            trimmed reads are mapped to the consensus of the previous iteration, or the provided reference genome in the
            first iteration. Afterwards, variants are called using <i>bcftools</i> (phase 1), followed by <i>Clair3</i>
            (phase 2). The workflow stops when the consensus sequence has converged (i.e., no changes in the last two
            iterations) or when the maximum number of iterations has passed (<b>{max_iter}</b>).
            """.format(max_iter=self._parameters['max_iter'].value))
        section.add_html_object(div)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        """
        section = HtmlReportSection('Consensus: iterative mapping')
        data_stats = pd.read_table(self._tool_inputs['TSV'][0].path)
        self.__add_section_iterative_mapping_stats(section, data_stats)
        self.__add_section_output_files(section, data_stats['dirname'].iloc[-1])
        self.__add_section_low_depth(section, data_stats)
        self.__add_section_additional_info(section)
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]
