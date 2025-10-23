import re
from pathlib import Path
from typing import Union

import pandas as pd

from camel.app.core.reports.htmlelement import HtmlElement
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class BTyperReporter(Tool):
    """
    This tool is used to generate HTML report sections based on the BTyper output.
    """

    TITLE = 'BTyper'
    URL_PUBMED = 'https://www.ncbi.nlm.nih.gov/pubmed/{id}'

    def __init__(self) -> None:
        """
        Initializes the tool.
        :return: None
        """
        super().__init__('BTyper Reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError('BTyper input (TSV) is required.')
        if 'btyper' not in self._input_informs:
            raise InvalidToolInputError('BTyper informs are required.')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the BTyper reporter tool
        :return: None
        """
        section = HtmlReportSection(BTyperReporter.TITLE, subtitle=self._input_informs['btyper']['_name'])

        # Add output tables
        formatted_tables = self.__parse_input_file()

        # Add species table
        section.add_header('Identification (ANI) analysis', level=3)
        species_table = formatted_tables['species']
        row_names = [HtmlElement('th', species_table.index.tolist()[k]) for k in range(len(species_table.index))]
        species_table_data = [[row_names[k], *species_table.values.tolist()[k]] for k in range(len(species_table.index))]
        section.add_table(species_table_data, ['', *species_table.columns], [('class', 'data')])

        # Add other tables
        self.__add_output_table(
            section, formatted_tables['virulence'].columns, formatted_tables['virulence'].values.tolist(), 'Virulence genes')
        self.__add_output_table(
            section, formatted_tables['typing'].columns, formatted_tables['typing'].values.tolist(), 'Typing')

        # Add link to TSV file
        relative_path = Path('btyper') / self._tool_inputs['TSV'][0].path.name
        section.add_link_to_file('Download complete results (TSV)', relative_path)
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path)

        # Tool output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]

    def __sanitize_table_lines(self, input_list: list[str]) -> list[str]:
        """
        Format nicely the entries in the HTML output table
        :param input_list: split line from the raw output table
        :return: List
        """
        for k in range(len(input_list)):
            if 'not performed' in input_list[k]:
                input_list[k] = 'Analysis not performed'

            # Remove empty brackets
            input_list[k] = re.sub(r'\(\)', '', input_list[k])

            # Add space before brackets
            input_list[k] = re.sub(r'(\w)\(', r'\1 (', input_list[k])

            #  Italics for gene names
            match_re = re.findall(r'[a-z]{3}[A-Z]', input_list[k])
            if match_re:
                input_list[k] = '{} ({})'.format(
                    input_list[k][:3], '; '.join([f'<i>{gene}</i>' for gene in match_re]))
        return input_list

    def _format_species_subtable(self, input_table: pd.DataFrame) -> pd.DataFrame:
        """
        Formats nicely the species subtable in two rows
        :param input_table: input species subtable
        :return: formatted species subtable
        """
        ani_regex = r'\d+\.\d+'
        species_regex = r'[a-z]+ ?[a-z]*\.?[a-z]*\.?'
        ani_species_match = re.findall(ani_regex, str(input_table['Species'].values[0]))
        if ani_species_match:
            species_match = re.findall(species_regex, str(input_table['Species'].values[0]))[0]
            ani_species_match = f'{float(ani_species_match[0]):.2f}'
        ani_subspecies_match = re.findall(ani_regex, str(input_table['Sub-species'].values[0]))
        if ani_subspecies_match:
            subspecies_match = re.findall(species_regex, str(input_table['Sub-species'].values[0]))[0]
            ani_subspecies_match = f'{float(ani_subspecies_match[0]):.2f}'
        table_dictionary = {'Species': ['n/a' if not ani_species_match else species_match,
                                        'n/a' if not ani_species_match else ani_species_match],
                            'Sub-species': ['n/a' if not ani_subspecies_match else subspecies_match,
                                            'n/a' if not ani_subspecies_match else ani_subspecies_match]}
        output_table = pd.DataFrame.from_dict(table_dictionary)
        output_table.index = ['Value', 'ANI (%)']
        return output_table

    def __parse_input_file(self) -> dict[str, pd.DataFrame]:
        """
        Parses the input file
        :return: Output tables
        """
        header_mapper = {'species(ANI)': 'Species', 'subspecies(ANI)': 'Sub-species',
                         'anthrax_toxin(genes)': 'Anthrax toxin genes',
                         'emetic_toxin_cereulide(genes)': 'Emetic toxin genes',
                         'diarrheal_toxin_Nhe(genes)': 'Diarrheal toxin <b>Nhe</b> genes',
                         'diarrheal_toxin_Hbl(genes)': 'Diarrheal toxin <b>Hbl</b> genes',
                         'diarrheal_toxin_CytK(top_hit)': 'Diarrheal toxin <b>CytK</b>',
                         'sphingomyelinase_Sph(gene)': 'Sphingomyelinase <b>Sph</b>',
                         'capsule_Cap(genes)': 'Capsule <b>Cap</b> genes',
                         'capsule_Has(genes)': 'Capsule <b>Has</b> genes',
                         'capsule_Bps(genes)': 'Capsule <b>Bps</b> genes', 'Bt(genes)': '<b>Bt</b> toxin genes',
                         'PubMLST_ST[clonal_complex](perfect_matches)': 'PubMLST ST',
                         'Adjusted_panC_Group(predicted_species)': 'Adjusted <i>panC</i> group',
                         'final_taxon_names': 'Final taxon name'}
        btyper_table = pd.read_table(self._tool_inputs['TSV'][0].path, sep='\t', usecols=list(header_mapper.keys()))
        btyper_table_values_sanitized = self.__sanitize_table_lines(btyper_table.values.tolist()[0])
        btyper_table = pd.DataFrame(btyper_table_values_sanitized).T
        btyper_table.columns = list(header_mapper.values())
        subtable_species = self._format_species_subtable(btyper_table[['Species', 'Sub-species']])
        subtable_virulence = btyper_table[
            ['Anthrax toxin genes', 'Emetic toxin genes', 'Diarrheal toxin <b>Nhe</b> genes',
             'Diarrheal toxin <b>Hbl</b> genes', 'Diarrheal toxin <b>CytK</b>',
             'Sphingomyelinase <b>Sph</b>', 'Capsule <b>Cap</b> genes', 'Capsule <b>Has</b> genes',
             'Capsule <b>Bps</b> genes', '<b>Bt</b> toxin genes']]
        subtable_typing = btyper_table[['PubMLST ST', 'Adjusted <i>panC</i> group', 'Final taxon name']]
        return {'species': subtable_species,
                'virulence': subtable_virulence,
                'typing': subtable_typing}

    def __generate_output_filename(self, prefix: str) -> str:
        """
        Generates the filename of the tabular output.
        :param prefix: prefix for the output file
        :return: Output filename
        """
        if 'sample_name' in self._parameters:
            return f"btyper-{self._parameters['sample_name'].value}_{prefix}.tsv"
        else:
            return f'btyper_{prefix}.tsv'

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
            section.add_header(f'{prefix} analysis', level=3)
            section.add_table(data, header, [('class', 'data')])
        else:
            section.add_paragraph('No results.')
