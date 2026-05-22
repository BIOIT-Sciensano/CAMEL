from pathlib import Path
from typing import Optional

import pandas as pd
from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltablecell import HtmlTableCell

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class AbriTAMRReporter(Tool):
    """
    Parses abriTAMR report output and generates an HTML output report.
    """

    TITLE = 'abriTAMR'

    def __init__(self) -> None:
        """
        Initializes the tool.
        """
        super().__init__('AbriTAMR Reporter', '0.1')
        self._section = None
        self._species = None

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if not all(key in self._tool_inputs for key in ('TXT_matches', 'TXT_partials')):
            raise InvalidToolInputError("AbriTAMR run output files must be provided")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the reporter.
        :rtype: None
        """
        self._section = HtmlReportSection(
            AbriTAMRReporter.TITLE,
            subtitle=self._input_informs['abritamr_run']['_name_full'],
        )
        self._species = self._input_informs['abritamr_run']['species']
        self.__add_summaries_tables()
        if self._species == 'Salmonella':
            self.__add_antibiogram()
        else:
            self._section.add_header('Antibiogram', 3)
            self._section.add_paragraph(f"Not available for species '{self._species}'")
        self.__add_output_table_link()
        self.__add_database_information()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_summaries_tables(self) -> None:
        """
        Adds summary tables for the antibiotic hit files not in the antibiogram.
        :return: None
        """
        self._section.add_header('Matches', 3)
        self._section.add_paragraph('>90% coverage & >90% identity')
        self.___add_summary_file_table(
            self._tool_inputs['TXT_matches'][0].path, 'matches'
        )
        self._section.add_paragraph(
            'Genes annotated with * indicate >90% coverage and identity between 90% and 100%. '
            'No further annotation indicates that the gene recovered exhibits '
            '100% coverage and 100% identity to a gene in the gene catalog.'
        )
        self._section.add_header('Partial matches', 3)
        self._section.add_paragraph('50-90% coverage & >90% identity')
        self.___add_summary_file_table(
            self._tool_inputs['TXT_partials'][0].path, 'partials'
        )
        self._section.add_horizontal_line()

    def ___add_summary_file_table(
        self, summary_file_path: Path, matching_type: str
    ) -> None:
        """
        Parses an AbriTAMR output summary file (TXT_partials or TXT_matches) and
        adds the table to the section.
        :param summary_file_path: path of TXT_partials or TXT_matches
        :param matching_type: str, matches or partials
        :return: None
        """
        header = ['Functional drug class', 'Genes']
        data = []

        # Read the summary file into a pandas DataFrame
        df = pd.read_table(summary_file_path)

        # Extract relevant data from DataFrame
        if len(df.columns) > 1:
            line1 = df.columns[1:]
            line2 = df.iloc[0, 1:].tolist()
            for functional_drug_class, gene in zip(line1, line2):
                data.append([functional_drug_class, gene])

        self._section.add_table(sorted(data), header, [('class', 'data')])
        relative_path_file = Path('abritamr', f"summary_{matching_type}.txt")
        self._section.add_file(summary_file_path, relative_path_file)
        self._section.add_link_to_file(
            f"Download {matching_type} (TSV)", relative_path_file
        )

    def __add_antibiogram(self) -> None:
        """
        Add the antibiogram to the report for Salmonella.
        :return: None
        """
        # Create usable dictionary from REPORT_abritamr for html_table
        antibiogram_dict = {}
        df_abritamr = pd.read_excel(
            self._tool_inputs['REPORT_abritamr'][0].path, engine='openpyxl', dtype=str
        )
        df_abritamr.fillna('-', inplace=True)  # replace all missing values by dashes
        for column in df_abritamr.columns[2:]:
            key_parts = column.replace(" - ", "_").split('_')
            value = df_abritamr.iloc[0][column]
            if key_parts[-1] in {'ResMech', 'Interpretation'}:
                antibiotic = '_'.join(
                    key_parts[:-1]
                )  # = join all except last value: counter intuitive
                antibiotic_property = key_parts[-1]

                # Initialize or update the sub dictionary in one step instead of checking with .get()
                antibiogram_dict.setdefault(antibiotic, {})[antibiotic_property] = value

        # Create html table
        header = ['Antibiotic', 'Resistance mechanisms detected', 'Interpretation']
        data = []
        for antibiotic, antibiotic_properties in antibiogram_dict.items():
            antibiotic_properties: dict[
                str, str
            ]  # add typing to silence PyCharm warnings
            color = self.___get_interpretation_color(
                antibiotic_properties['Interpretation']
            )
            row = [
                antibiotic,
                antibiotic_properties['ResMech'],
                antibiotic_properties['Interpretation'],
            ]
            row = [HtmlTableCell(x, color) for x in row]
            data.append(row)
        # start writing in the report the table and the headers
        self._section.add_header('Antibiogram for <i>Salmonella</i>', 3)
        self._section.add_table(data, header, [('class', 'data')])

    @staticmethod
    def ___get_interpretation_color(interpretation: str) -> Optional[str]:
        """
        For a given interpretation returns either green, red or None
        :param interpretation: Resistant or Susceptible or any other string
        :return: str ('green' or 'red') or None
        """
        if interpretation == 'Susceptible':
            color = 'green'
        elif interpretation == 'Resistant':
            color = 'red'
        else:
            color = 'grey'
        return color

    def __add_output_table_link(self) -> None:
        """
        Add link to the output table (tsv) for this assay.
        :return: None
        """
        if self._species == 'Salmonella':
            relative_path = Path('abritamr', 'summary_out.xlsx')
            self._section.add_file(
                self._tool_inputs['REPORT_abritamr'][0].path, relative_path
            )
            self._section.add_link_to_file("Download (xlsx)", relative_path)

    def __add_database_information(self) -> None:
        """
        Adds the database information to the report.
        :return: None
        """
        self._section.add_horizontal_line()
        self._section.add_header('Database info', level=4)
        self._section.add_paragraph(
            f'Last updated: {self._input_informs["abritamr_run"].get("last_update_date", "n/a")}'
        )
