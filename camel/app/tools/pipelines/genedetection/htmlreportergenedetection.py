from pathlib import Path

from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlelement import HtmlElement
from camelcore.app.reports.htmlexpandablediv import HtmlExpandableDiv
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltablecell import HtmlTableCell
from camelcore.app.utils import fileutils

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool
from camel.app.loggers import logger
from camel.app.toolkits.export.tsvexporter import TsvExporter
from camel.app.toolkits.genedetection.genedetectionhitbase import GeneDetectionHitBase
from camel.app.toolkits.genedetection.genedetectionmixedhit import GeneDetectionMixedHit


class HtmlReporterGeneDetection(Tool):
    """
    Tool that creates HTML reports for the gene detection pipeline.
    """

    COLS_MIXED = {
        'blast': [{
            'name': '% identity',
            'func': lambda x: f'{x.blast_stats.subject_normalized_identity:.2f}'
        }, {
            'name': '% covered',
            'func': lambda x: f'{x.blast_stats.subject_coverage:.2f}'
        }],
        'kma': [
            {
                'name': '% identity',
                'func': lambda x: f'{x.percent_identity:.2f}'
            },
            {
                'name': '% covered',
                'func': lambda x: f'{x.subject_coverage:.2f}'
            },
            {
                'name': 'Depth',
                'func': lambda x: f'{x.depth:,}'
            },
        ],
    }

    def __init__(self) -> None:
        """
        Initialize this tool.
        :return: None
        """
        super().__init__('Gene Detection: Report', '0.1')
        self._sub_folder = None
        self._report_section = None

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__initialize_report()

        # Parameter section
        self.__add_parameter_table(self._input_informs['detection'])

        # Detected hits
        if len(self._tool_inputs['VAL_Hits']) == 0:
            self._report_section.add_paragraph('No hits found.')
        else:
            self.__add_output_table([h.value for h in self._tool_inputs['VAL_Hits']])

        # For mixed gene detection, show the table with additional details
        if self.get_param_value('detection_method') == 'mixed':
            self.__add_table_mixed_long([io.value for io in self._tool_inputs['VAL_Hits']])

        # Additional information
        if 'extra_info' in self._input_informs['detection']:
            self.__add_additional_info_section_aggregated_hits(self._input_informs['detection']['extra_info'])

        # Status messages
        self.__add_messages()

        # Database information
        self.__add_database_information()

        # Create output
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section)]

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        if 'db_info' not in self._input_informs:
            raise InvalidToolInputError("No database info found")
        if 'VAL_Hits' not in self._tool_inputs:
            logger.warning("No BLAST hits found")
        if ('VAL_Hits' in self._tool_inputs) and (len(self._tool_inputs['VAL_Hits']) > 0) and \
                ('TSV' not in self._tool_inputs):
            raise InvalidToolInputError("TSV input is required when hits were detected.")
        super()._check_input()

    def __initialize_report(self) -> None:
        """
        Initializes the HTML report.
        :return: None
        """
        db_name = self._input_informs['db_info']['name']
        self._report_section = HtmlReportSection(self._input_informs['db_info']['title'], 3)
        self._sub_folder = Path('gene_detection') / fileutils.make_valid(db_name)

    def __add_parameter_table(self, informs_detection: dict[str, str]) -> None:
        """
        Adds a tables with the parameters used for the detection.
        :param informs_detection: Informs from the detection
        :return: None
        """
        self._report_section.add_table([
            [f'{key}:', value] for key, value in sorted(informs_detection.items()) if not key.startswith('_')
            if not key.startswith('_')
            and key != 'extra_info'
        ], None, [('class', 'information')])

    def __add_output_table(self, hits: list[GeneDetectionHitBase]) -> None:
        """
        Adds the output table.
        :param hits: Detected hits
        :return: None
        """
        table_data = [hit.to_html_row(self._report_section, self._sub_folder) for hit in sorted(
            hits, key=lambda x: x.locus)]
        if self.get_param_value('hidden') is True:
            div = HtmlExpandableDiv('table-{}'.format(
                self._input_informs['db_info']['name']), f'hits ({len(hits):,})')
            div.add_table(table_data, hits[0].html_column_names, [('class', 'data')])
            self._report_section.add_html_object(div)
        else:
            self._report_section.add_table(table_data, hits[0].html_column_names, [('class', 'data')])
        relative_path = self._sub_folder / Path(self._tool_inputs['TSV'][0].path).name
        self._report_section.add_file(self._tool_inputs['TSV'][0].path, relative_path)
        self._report_section.add_link_to_file("Download (TSV)", relative_path)

    @staticmethod
    def __first_locus_name(hit_mixed: GeneDetectionMixedHit) -> str:
        """
        Returns the alphabetically first locus name detected for this hit's cluster, used to sort clusters.
        :param hit_mixed: Mixed hit
        :return: Alphabetically first locus name, or '' if no loci were detected
        """
        locus_names = [(locus_hits['kma'] or locus_hits['blast']).locus for locus_hits in hit_mixed.hits_by_locus.values()]
        return min(locus_names, default='')

    def __add_table_mixed_long(self, hits: list[GeneDetectionMixedHit]) -> None:
        """
        Adds the complete list of detected hits as a "long" HTML table, with one row per locus. A locus detected
        by both KMA and Blast is rendered on a single row with both tools' stats; a locus detected by only one
        tool is rendered with a "Not detected by <other tool>" placeholder. The cluster cell is merged (rowspan)
        across all locus rows belonging to the same cluster. Rows are sorted alphabetically by locus name, both
        within a cluster and across clusters (a cluster is ordered by its alphabetically first locus).
        :param hits: Mixed hits to render, one per cluster
        :return: None
        """
        self._report_section.add_header('All detected hits (complete list)', level=3)
        self._report_section.add_paragraph(
            'This table lists every locus detected by KMA and/or BLAST, along with the detailed statistics from each '
            'method.')

        ##########
        # Header #
        ##########
        tsv_header = [
            'Cluster', 'Locus', 'Accession', 'Length',
            *(f"BLAST {info['name']}" for info in HtmlReporterGeneDetection.COLS_MIXED['blast']),
            *(f"KMA {info['name']}" for info in HtmlReporterGeneDetection.COLS_MIXED['kma']),
        ]
        header = [
            [
                HtmlElement('th', 'Cluster', [('rowspan', 2)]),
                HtmlElement('th', 'Locus', [('rowspan', 2)]),
                HtmlElement('th', 'Accession', [('rowspan', 2)]),
                HtmlElement('th', 'Length', [('rowspan', 2)]),
                HtmlElement('th', 'BLAST', [('colspan', len(HtmlReporterGeneDetection.COLS_MIXED['blast']))]),
                HtmlElement('th', 'KMA', [('colspan', len(HtmlReporterGeneDetection.COLS_MIXED['kma']))]),
            ],
            [
                *(HtmlElement('th', info['name']) for info in HtmlReporterGeneDetection.COLS_MIXED['blast']),
                *(HtmlElement('th', info['name']) for info in HtmlReporterGeneDetection.COLS_MIXED['kma']),
            ]
        ]

        ##############
        # Table rows #
        ##############
        table_data = []
        tsv_data = []
        for hit_mixed in sorted(hits, key=self.__first_locus_name):
            rows_for_cluster = []
            locus_hits_sorted = sorted(
                hit_mixed.hits_by_locus.values(), key=lambda locus_hits: (locus_hits['kma'] or locus_hits['blast']).locus)
            for locus_hits in locus_hits_sorted:
                hit_kma, hit_blast = locus_hits['kma'], locus_hits['blast']
                locus_hit = hit_kma or hit_blast
                length = hit_kma.length if hit_kma is not None else hit_blast.blast_stats.subject_length
                row = [
                    HtmlElement('td', locus_hit.locus),
                    locus_hit._get_accession_cell(colored=False),
                    HtmlElement('td', length),
                ]
                if hit_blast is not None:
                    row += [
                        HtmlTableCell(info['func'](hit_blast), hit_blast.color) for info in HtmlReporterGeneDetection.COLS_MIXED['blast']
                    ]
                else:
                    row.append(HtmlElement(
                        'td', 'Not detected by BLAST', [('colspan', len(HtmlReporterGeneDetection.COLS_MIXED['blast']))]))
                if hit_kma is not None:
                    row += [
                        HtmlTableCell(info['func'](hit_kma), hit_kma.color) for info in HtmlReporterGeneDetection.COLS_MIXED['kma']
                    ]
                else:
                    row.append(HtmlElement(
                        'td', 'Not detected by KMA', [('colspan', len(HtmlReporterGeneDetection.COLS_MIXED['kma']))]))
                rows_for_cluster.append(row)
                tsv_data.append([
                    hit_mixed.cluster,
                    locus_hit.locus,
                    locus_hit.accession if locus_hit.accession is not None else '-',
                    length,
                    *(info['func'](hit_blast) if hit_blast is not None else '-'
                        for info in HtmlReporterGeneDetection.COLS_MIXED['blast']),
                    *(info['func'](hit_kma) if hit_kma is not None else '-'
                        for info in HtmlReporterGeneDetection.COLS_MIXED['kma']),
                ])
            if not rows_for_cluster:
                continue

            # Merge the cluster cell across all rows of this cluster, then add the remaining rows as-is
            table_data.append(
                [HtmlElement('td', hit_mixed.cluster, [('rowspan', len(rows_for_cluster))])] + rows_for_cluster[0]
            )
            table_data.extend(rows_for_cluster[1:])

        # Add the table, in an expandable div
        div = HtmlExpandableDiv('table_complete_agg_hits', f'detected hits ({len(tsv_data):,})')
        div.add_table(header + table_data, None, [('class', 'data')])
        self._report_section.add_html_object(div)

        # Create TSV file and add download link
        output_tsv = Path(self._folder) / 'agg-hits-summary.tsv'
        TsvExporter.export(tsv_data, tsv_header, str(output_tsv))
        relative_path = self._sub_folder / output_tsv.name
        self._report_section.add_file(output_tsv, relative_path)
        self._report_section.add_link_to_file("Download (TSV)", relative_path)

    def __add_additional_info_section_aggregated_hits(self, extra_info: str) -> None:
        """
        Adds the complete aggregated hit list section to the report.
        :param extra_info: paragraph to add
        :return: None
        """
        self._report_section.add_header('Additional information', level=4)
        self._report_section.add_paragraph(extra_info)

    def __add_messages(self) -> None:
        """
        Adds status messages to the report: a notice when the detection method is forced, a warning when run on
        simulated reads, and a custom message (if configured).
        :return: None
        """
        if 'forced_detection_method' in self._parameters:
            self._report_section.add_alert(
                f"Detection for this DB is always done using '{self._parameters['forced_detection_method'].value}', "
                f"regardless of pipeline setting.", 'info')

        if self.get_param_value('pseudo_reads') is True:
            self._report_section.add_warning_message("The tool is executed on simulated reads.")

        if 'message' in self._parameters:
            self._report_section.add_alert(
                self._parameters['message'].value, self._parameters['message_category'].value)

    def __add_database_information(self) -> None:
        """
        Adds the database information to the report.
        :return: None
        """
        self._report_section.add_header('Database info', level=4)
        self._report_section.add_table([
            ('Last database update', self._input_informs['db_info']['last_updated']),
            ('Last database change', self._input_informs['db_info'].get('last_change', 'n/a'))
        ], ['Field', 'Value'], [('class', 'data')])
