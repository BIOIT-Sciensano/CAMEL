from pathlib import Path
from typing import Union

import pandas as pd

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class StrainGSTReporter(Tool):
    """
    This tool is used to generate HTML report sections based on the StrainGST output.
    """

    TITLE = 'StrainGST'
    URL_PUBMED = 'https://www.ncbi.nlm.nih.gov/pubmed/{id}'

    def __init__(self) -> None:
        """
        Initializes the tool.
        """
        super().__init__('StrainGST Reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError('StrainGST strains input (TSV) is required.')
        if 'straingst' not in self._input_informs:
            raise InvalidToolInputError('StrainGST informs are required.')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the StrainGST reporter tool
        :return: None
        """
        section = HtmlReportSection(StrainGSTReporter.TITLE, subtitle=self._input_informs['straingst']['_name'])
        suffix_read_type = self._parameters["suffix"].value.upper()

        # Add output tables
        output_table = self.__parse_input_file()
        self.__add_output_table(
            section, list(output_table.columns), output_table.values.tolist(),
            f'StrainGST strain identification - {suffix_read_type}')

        # Add link to TSV file
        relative_path = Path('straingst') / self._tool_inputs['TSV'][0].path.name
        section.add_link_to_file('Download complete results (TSV)', relative_path)
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path)

        # Tool output
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]

    def __parse_input_file(self) -> pd.DataFrame:
        """
        Parses the input file.
        :return: Dataframe with relevant columns
        """
        to_return = pd.read_table(
            self._tool_inputs['TSV'][0].path, header=0,
            usecols=[1, 5, 9, 11, 14],
            names=['Strain', 'Coverage', 'Evenness', 'Relative abundance', 'Score'])
        to_return['Coverage'] = pd.Series([
            f"{float(val) * 100:.2f}%" for val in to_return['Coverage']], index=to_return.index)
        to_return['Evenness'] = pd.Series(
            [f"{float(val) * 100:.2f}%" for val in to_return['Evenness']], index=to_return.index)
        to_return['Relative abundance'] = pd.Series(
            [f"{float(val):.2f}%" for val in to_return['Relative abundance']], index=to_return.index)
        return to_return

    def __generate_output_filename(self, prefix: str) -> str:
        """
        Generates the filename of the tabular output.
        :return: Output filename
        """
        if 'sample_name' in self._parameters:
            return f"straingst-{self._parameters['sample_name'].value}_{prefix}.tsv"
        else:
            return f'straingst_{prefix}.tsv'

    def __add_output_table(
            self, section: HtmlReportSection, header: list[str],
            data: list[list[Union[str, HtmlTableCell]]], prefix: str) -> None:
        """
        Adds an output table to the HTML report.
        :param section: Report section
        :param header: Output table header
        :param data: Output table data
        :param prefix: Prefix for the table to add
        :return: None
        """
        if len(data) > 0:
            section.add_header(f'{prefix}', level=3)
            section.add_table(data, header, [('class', 'data')])
        else:
            section.add_paragraph('No results.')
