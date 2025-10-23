from pathlib import Path
from typing import Any, Optional

import pandas as pd

from camel.app.core.reports.htmlelement import HtmlElement
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.core.tool import Tool


class HtmlReporterContamination(Tool):
    """
    HTML reporter for the contamination check step.
    """

    TITLE = 'Contamination check'

    ABBREVIATIONS = {
        'U': 'Unclassified',
        'D': 'Domain',
        'K': 'Kingdom',
        'P': 'Phylum',
        'C': 'Class',
        'O': 'Order',
        'F': 'Family',
        'G': 'Genus',
        'S': 'Species'
    }

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('HTML Reporter: Contamination', '0.1')
        self._report_section = None

    @property
    def _dir_out(self) -> Path:
        """
        Returns the relative path to the output directory.
        :return: Relative path to output directory
        """
        if 'suffix' not in self._parameters:
            return Path('contamination_check')
        return Path('contamination_check', str(self._parameters['suffix'].value))

    @property
    def _title(self) -> str:
        """
        Returns the title for the output report.
        :return: Report title
        """
        if 'suffix_title' not in self._parameters:
            return 'Contamination check'
        return f"Contamination  check - {self._parameters['suffix_title'].value}"

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._report_section = HtmlReportSection(
            self._title, subtitle=self._input_informs['kraken2']['_name'])
        self.__add_database_info()
        self.__add_filtering_info(self._input_informs['species'])
        self.__add_species_table()
        logger.info(self._input_informs)
        self.__add_detailed_table(self._tool_inputs['TSV'][0].path)
        if self._input_informs['species']['normalize_for_len']:
            self._report_section.add_alert(
                '(*) These percentages have been normalised to take into account the length of the reads / contigs.'
                ' This normalisation is not performed for the Krona report.', 'info')
        self.__add_krona_report()
        self.__add_warnings()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section)]

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'HTML_Krona' not in self._tool_inputs:
            raise InvalidToolInputError("Krona report input is required.")
        if 'species' not in self._input_informs:
            raise InvalidToolInputError("Species input is required")
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError("Kraken report input (TSV) is required")
        super()._check_input()

    def __add_database_info(self) -> None:
        """
        Adds the database information to the report.
        :return: None
        """
        db_informs = self._input_informs['kraken2'].get('database')
        self._report_section.add_table([
            ['Database:', db_informs['name'] if db_informs is not None else 'NA'],
            ['Last update:', db_informs['last_update'] if db_informs is not None else 'NA']
        ], None, [('class', 'information')])

    def __add_filtering_info(self, informs: dict[str, Any]) -> None:
        """
        Adds the filtering parameters.
        :param informs: species informs
        :return: None
        """
        self._report_section.add_table([
            ['Warning threshold:', f"{informs['threshold_warn']}%"],
            ['Fail threshold:', f"{informs['threshold_fail']}%"],
            ['Level:', '{}'.format('Species' if 'S' in informs['level_of_depth'] else 'Genus')]
        ], None, [('class', 'information')])

    def __add_species_table(self) -> None:
        """
        Adds a table containing the detected species and corresponding percentages.
        :return: None
        """
        header = ['Species' if 'S' in self._input_informs['species']['level_of_depth'] else 'Genus', 'Percentage']
        expected_name, expected_perc = self._input_informs['species']['expected']
        expected_perc = f'{float(expected_perc):.2f}'
        table_data = [
            [HtmlElement('th', 'Expected', [('colspan', 2)])],
            [HtmlTableCell(f'<i>{expected_name}</i>', 'green'), HtmlTableCell(expected_perc, 'green')],
        ]
        # Allowed species
        if len(self._input_informs['species']['allowed']) > 0:
            table_data.append([HtmlElement('th', 'Allowed', attributes=[('colspan', 2)])])
            for species_name, percentage in self._input_informs['species']['allowed']:
                table_data.append([
                    HtmlTableCell(f'<i>{species_name}</i>', 'yellow'),
                    HtmlTableCell(f'{percentage:.2f}', 'yellow')
                ])

        # Failed contaminants
        table_data.append([HtmlElement('th', 'Contaminants', [('colspan', 2)])])
        if (len(self._input_informs['species']['contaminants_fail']) +
                len(self._input_informs['species']['contaminants_warn']) == 0):
            table_data.append([HtmlTableCell('None found', attributes=[('colspan', 2)])])
        for species_name, percentage in self._input_informs['species']['contaminants_fail']:
            table_data.append([
                HtmlTableCell(f'<i>{species_name}</i>', 'red'),
                HtmlTableCell(f'{percentage:.2f}', 'red')])
        # Warning contaminants
        for species_name, percentage in self._input_informs['species']['contaminants_warn']:
            table_data.append([
                HtmlTableCell(f'<i>{species_name}</i>', 'orange'),
                HtmlTableCell(f'{percentage:.2f}', 'orange')])

        self._report_section.add_table(table_data, header, [('class', 'data')])

    def __add_detailed_table(self, kraken_report_path: Path) -> None:
        """
        Adds a table with the detailed information from the tabular KRAKEN output.
        :return: None
        """
        data_report = pd.read_table(kraken_report_path)
        perc_is_normalized = 'perc_orig' in data_report.columns
        table_data = []
        for row in data_report.to_dict('records'):
            if row['perc'] < 1.0:
                continue
            name = f'<i>{row["name"].strip()}</i>' if row['rank'] in ('G', 'S') else row['name'].strip()
            table_data.append([
                f"{row['perc']:.2f}",
                HtmlReporterContamination.ABBREVIATIONS.get(row['rank'], '-'),
                name])
            if perc_is_normalized:
                table_data[-1].insert(1, row['perc_orig'])
        header = ['Percentage', 'Level', 'Name'] if not perc_is_normalized else ['Perc. (*)', 'Perc. orig.', 'Level', 'Name']
        self._report_section.add_table(table_data, header, [('class', 'data')])

    def __add_krona_report(self) -> None:
        """
        Adds a download link to the Krona report.
        :return: None
        """
        relative_path = self._dir_out / 'krona_report.html'
        self._report_section.add_file(self._tool_inputs['HTML_Krona'][0].path, relative_path)
        self._report_section.add_link_to_file('Krona Report', relative_path)
        self._report_section.add_file(self._tool_inputs['TSV'][0].path, self._dir_out / 'kraken2_report.tsv')

    def __add_warnings(self) -> None:
        """
        Adds warnings to the report when there is possible contamination.
        :return: None
        """
        if len(self._input_informs['species']['contaminants_warn']) > 0:
            self._report_section.add_warning_message("The sample is possibly contaminated with: {}".format(
                ", ".join(f'<i>{species_name}</i>' for species_name, _ in
                          self._input_informs['species']['contaminants_warn'])))
        if len(self._input_informs['species']['contaminants_fail']) > 0:
            self._report_section.add_error_message("There is strong evidence of a contamination with: {}".format(
                ", ".join(f'<i>{species_name}</i>' for species_name, _ in
                          self._input_informs['species']['contaminants_fail'])))

    @staticmethod
    def generate_empty_section(suffix_title: Optional[str] = None) -> HtmlReportSection:
        """
        Returns a report used when this analysis is disabled.
        :param suffix_title: (Optional) suffix for the title
        :return: Report section
        """
        section = HtmlReportSection(
            HtmlReporterContamination.TITLE if suffix_title is None else
            f'{HtmlReporterContamination.TITLE} - {suffix_title}')
        section.add_paragraph('Analysis disabled')
        return section
