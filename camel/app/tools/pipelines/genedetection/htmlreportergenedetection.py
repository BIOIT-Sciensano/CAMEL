from pathlib import Path

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.reports.htmlexpandablediv import HtmlExpandableDiv
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.tool import Tool
from camel.app.core.utils import fileutils
from camel.app.loggers import logger
from camel.app.toolkits.genedetection.genedetectionhitbase import GeneDetectionHitBase
from camel.app.tools.pipelines.genedetection.mixedhitaggregating import AggregatedHit


class HtmlReporterGeneDetection(Tool):
    """
    Tool that creates HTML reports for the gene detection pipeline.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
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
        self.__add_parameter_table(self._input_informs['detection'])
        if len(self._tool_inputs['VAL_Hits']) == 0:
            self._report_section.add_paragraph('No hits found.')
        else:
            self.__add_output_table([h.value for h in self._tool_inputs['VAL_Hits']])
            if self._tool_inputs['VAL_agg_hits']:
                self.__add_additional_info_section_aggregated_hits()

        self.__add_database_information()

        # Add a warning when the detection method is different from the general detection
        if 'forced_detection_method' in self._parameters:
            self._report_section.add_alert(
                f"Detection for this DB is always done using '{self._parameters['forced_detection_method'].value}', "
                f"regardless of pipeline setting.", 'info')

        if 'pseudo_reads' in self._parameters:
            self._report_section.add_warning_message("The tool is executed on simulated reads.")

        # Add a custom message (if specified in the parameters)
        if 'message' in self._parameters:
            self._report_section.add_alert(
                self._parameters['message'].value, self._parameters['message_category'].value)

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
            logger.warning("No blast hits found")
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
        ], None, [('class', 'information')])

    def __add_output_table(self, hits: list[GeneDetectionHitBase]) -> None:
        """
        Adds the output table.
        :param hits: Detected hits
        :return: None
        """
        table_data = [hit.to_html_row(self._report_section, self._sub_folder) for hit in sorted(
            hits, key=lambda x: x.locus)]
        if 'hidden' in self._parameters:
            div = HtmlExpandableDiv('table-{}'.format(
                self._input_informs['db_info']['name']), f'hits ({len(hits):,})')
            div.add_table(table_data, hits[0].html_column_names, [('class', 'data')])
            self._report_section.add_html_object(div)
        else:
            self._report_section.add_table(table_data, hits[0].html_column_names, [('class', 'data')])
        relative_path = self._sub_folder / Path(self._tool_inputs['TSV'][0].path).name
        self._report_section.add_file(self._tool_inputs['TSV'][0].path, relative_path)
        self._report_section.add_link_to_file("Download (TSV)", relative_path)

    def __add_aggregated_table(self, hits: list[AggregatedHit]) -> None:
        """
        Adds the aggregated table to the report.
        This table is added only if the "mixed" gene detection approach is selected
        and includes all genes detected by either method.
        :param hits: Detected aggregated hits
        :return: None
        """
        self._report_section.add_table([hit.to_html_list() for hit in sorted(hits, key=lambda x: x.locus)],
                                       column_names=hits[0].header()[:-2],
                                       table_attributes=[('class', 'data')])
        relative_path = self._sub_folder / Path(self._tool_inputs['TSV_agg_hits'][0].path).name
        self._report_section.add_file(self._tool_inputs['TSV_agg_hits'][0].path, relative_path)
        self._report_section.add_link_to_file("Download (TSV)", relative_path)

    def __add_additional_info_section_aggregated_hits(self) -> None:
        """
        Adds the complete aggregated hit list section to the report.
        :return: None
        """
        self._report_section.add_header('Additional information', level=4)
        self._report_section.add_paragraph(
            'In the case of the mixed gene detection strategy, both <b>read-based (KMA)</b> and <b>assembly-based (Blast)</b> '
            'gene detection methods are used. '
            'Since Blast can report hits on multiple closely related alleles, and KMA only single best hits, '
            'the table contains matching alleles detected with both KMA and Blast. '
            'The details about all detected genes (i.e., all alleles) are presented in the following table.'
        )
        self._report_section.add_header('Aggregated hits: complete list.', level=3)
        self._report_section.add_paragraph(
            'This table shows the details about all genes detected by either strategies.'
        )
        self.__add_aggregated_table([h.value for h in self._tool_inputs['VAL_agg_hits']])

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
