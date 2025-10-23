from pathlib import Path

import pandas as pd

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class AMRFinderReporter(Tool):
    """
    Reporting class for the AMRFinder tool.
    """

    URL_ACCESSION_BASE = 'https://www.ncbi.nlm.nih.gov/protein/'
    OUTPUT_COLS_OVERVIEW = [
        'Element symbol', 'Type', 'Subtype', 'Class', 'Subclass', 'Method',
        'Closest reference accession']
    OUTPUT_COLS_ALN = [
        'Element symbol', 'Contig id', 'Start', 'Stop', 'Strand', 'Target length', 'Alignment length',
        'Reference sequence length', '% Coverage of reference', '% Identity to reference']
    OUTPUT_COL_ACCESSION = 'Closest reference accession'

    def __init__(self) -> None:
        """
        Initializes the tool.
        :return: None
        """
        super().__init__('AMRFinder reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError("TSV input is required")
        if 'amrfinder' not in self._input_informs:
            raise InvalidToolInputError("AMRFinder informs are required")
        super()._check_input()

    def __add_stats_table(self, section: HtmlReportSection, df_in: pd.DataFrame, columns: list[str]) -> None:
        """
        Formats the column data.
        :param section: Report section
        :param df_in: Input dataframe
        :param columns: Selected columns
        :return: None
        """
        # Reformat data
        rows_out = []
        for index, row in df_in.iterrows():
            row_current = []
            for col in columns:
                if col.startswith('%'):
                    # XX.XX% notation for percentages
                    row_current.append(HtmlTableCell(f'{row[col]:.2f}', color=row['row_color']))
                elif col == AMRFinderReporter.OUTPUT_COL_ACCESSION:
                    # Links for accession numbers
                    row_current.append(HtmlTableCell(
                        row[col], link=f'{AMRFinderReporter.URL_ACCESSION_BASE}{row[col]}', color=row['row_color']))
                else:
                    row_current.append(HtmlTableCell(row[col], color=row['row_color']))
            rows_out.append(row_current)

        # Add table
        section.add_table(rows_out, columns, [('class', 'data')])

    @staticmethod
    def __determine_row_color(row_data: pd.Series) -> str:
        """
        Determines the row color.
        :param row_data: Row data
        :return: Color
        """
        if all(row_data[x] == 100.0 for x in ('% Identity to reference', '% Coverage of reference')):
            return 'green'
        elif row_data['% Coverage of reference'] == 100.0:
            return 'lightgreen'
        else:
            return 'grey'

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('AMRFinder', subtitle=self._input_informs['amrfinder']['_name_full'])

        # Parse input data
        data_out = pd.read_table(self._tool_inputs['TSV'][0].path)

        # Overview table
        if len(data_out) > 0:
            data_out['row_color'] = data_out.apply(lambda x: AMRFinderReporter.__determine_row_color(x), axis=1)
            section.add_header('Overview', level=3)
            data_out['Class'] = data_out['Class'].apply(lambda x: ', '.join([p.title() for p in x.split('/')]))
            data_out['Subclass'] = data_out['Subclass'].apply(lambda x: ', '.join([p.title() for p in x.split('/')]))
            self.__add_stats_table(section, data_out, AMRFinderReporter.OUTPUT_COLS_OVERVIEW)

            # Alignment stats table
            section.add_header('Alignment', level=3)
            self.__add_stats_table(section, data_out, AMRFinderReporter.OUTPUT_COLS_ALN)
        else:
            section.add_paragraph('No hits found')

        # Output and info section
        section.add_header('Output and extra information', level=3)
        relative_path = Path('amrfinder') / self._tool_inputs['TSV'][0].path.name
        section.add_link_to_file('Download complete results (TSV)', relative_path)
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path)
        section.add_paragraph(f"Database version: {self._input_informs['amrfinder']['db_version']}")

        # Set the tool outputs
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
