from pathlib import Path

from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.genedetection.genedetectionhitbase import GeneDetectionHitBase
from camel.app.components.html.htmlexpandablediv import HtmlExpandableDiv
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


class HtmlReporterGeneDetection(Tool):
    """
    Tool that creates HTML reports for the gene detection pipeline.
    """

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
        self.__add_parameter_table(self._input_informs['detection'])
        if len(self._tool_inputs['VAL_Hits']) == 0:
            self._report_section.add_paragraph('No hits found.')
        else:
            self.__add_output_table([h.value for h in self._tool_inputs['VAL_Hits']])
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
        self._sub_folder = Path('gene_detection') / FileSystemHelper.make_valid(db_name)

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
