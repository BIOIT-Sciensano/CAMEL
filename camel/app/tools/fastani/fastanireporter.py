from pathlib import Path
from typing import Union

import pandas as pd

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class FastANIReporter(Tool):
    """
    This tool is used to generate HTML report sections based on the FastANI output.
    """

    TITLE = 'FastANI'
    URL_PUBMED = 'https://www.ncbi.nlm.nih.gov/pubmed/{id}'

    def __init__(self) -> None:
        """
        Initializes the tool.
        """
        super().__init__('FastANI Reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError('FastANI input (TSV) is required.')
        if 'fastani' not in self._input_informs:
            raise InvalidToolInputError('FastANI informs are required.')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the FastANI reporter tool
        :return: None
        """
        section = HtmlReportSection(FastANIReporter.TITLE, subtitle=self._input_informs['fastani']['_name'])

        # Add output tables
        output_table = self.__parse_input_file()
        self.__add_output_table(section, output_table.columns, output_table.values.tolist())

        # Add link to TSV file
        relative_path = Path('fastani') / self._tool_inputs['TSV'][0].path.name
        section.add_link_to_file('Download complete results (TSV)', relative_path)
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path)

        # Tool output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]

    def __parse_input_file(self) -> pd.DataFrame:
        """
        Parses the FastANI input file.
        :return: formatted dataframe
        """
        with open(self._tool_inputs['TSV'][0].path) as handle:
            header = ['Query', 'ANI', 'Orthologous matches', 'Total seq fragments']
            fastani_table = pd.read_table(handle, sep='\t', header=0, names=header)
            fastani_table['ANI (%)'] = fastani_table['ANI'].map("{:.2f}".format)
            fastani_table['Query'] = fastani_table.apply(lambda row: self.__format_output_table_line(row['Query']),
                                                         axis=1)
            return fastani_table.head(10)

    def __format_output_table_line(self, query: str) -> str:
        """
        Formats a line in the output table to have only filenames and 2 significant digits on ANI.
        :query: query column value from the fastani table.
        :return: formatted genome name
        """
        if 'species' in self._parameters and self._parameters['species'].value == 'subtilis':
            # return a clean species name if bacillus is used (well formatted table)
            species = Path(str(query)).parent.parent.parent.name
            subspecies = Path(str(query)).parent.parent.name
            strain = Path(str(query)).parent.name
            query_out = f'{species}_{subspecies}_str_{strain}'
        else:
            # Remove directories from filenames
            query_out = Path(str(query)).name
        return query_out

    def __generate_output_filename(self) -> str:
        """
        Generates the filename of the tabular output.
        :return: Output filename
        """
        if 'sample_name' in self._parameters:
            return f"fastani-{self._parameters['sample_name'].value}.tsv"
        else:
            return 'fastani.tsv'

    def __add_output_table(
            self, section: HtmlReportSection, header: list[str],
            data: list[list[Union[str, HtmlTableCell]]]) -> None:
        """
        Adds an output table to the HTML report.
        :param section: Report section
        :param header: Output table header
        :param data: Output table data
        :return: None
        """
        if len(data) > 0:
            section.add_header('FastANI results', level=3)
            section.add_table(data, header, [('class', 'data')])
        else:
            section.add_paragraph('No results.')
