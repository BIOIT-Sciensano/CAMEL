import itertools
import json
import re
from pathlib import Path
from typing import Any, Optional

from camel.app.core.reports.htmlelement import HtmlElement
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.tools.pipelines.mycobacterium.amr import amrtypedetermination
from camel.app.core.tool import Tool
from camel.app.toolkits.export.tsvexporter import TsvExporter
from camel.app.toolkits.mycobacterium import amrutils


class AMRReporter(Tool):
    """
    Creates reports for the AMR detection for M. tuberculosis.
    """

    TITLE = 'AMR detection'

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Mycobacterium: AMR reporter', '0.1')
        self._sub_folder = 'amr'
        self._section = HtmlReportSection(AMRReporter.TITLE)
        self._actg_counts_by_pos = None
        self._mutations = None
        self._ab_info_by_category = {}

    def _check_input(self) -> None:
        """
        Checks if the required input is present.
        :return: None
        """
        if 'JSON' not in self._tool_inputs:
            raise InvalidToolInputError("AMR association info is required ('JSON')")
        if 'JSON_counts' not in self._tool_inputs:
            raise InvalidToolInputError("Nucleotide counts are required ('JSON_counts')")
        if 'JSON_pheno' not in self._tool_inputs:
            raise InvalidToolInputError("Predicted phenotypes are required ('JSON_pheno')")
        if 'JSON_amr_type' not in self._tool_inputs:
            raise InvalidToolInputError("AMR type information is required ('JSON_amr_type')")
        if 'PNG' not in self._tool_inputs:
            raise InvalidToolInputError("Circos visualization is required ('PNG')")
        if 'screen' not in self._input_informs:
            raise InvalidToolInputError("Mutation screening informs are required")

    def __parse_input_files(self) -> None:
        """
        Parses the input files.
        :return: None
        """
        with open(self._tool_inputs['JSON_counts'][0].path) as handle:
            self._actg_counts_by_pos = {int(pos): counts for pos, counts in json.load(handle).items()}

        with open(self._tool_inputs['JSON'][0].path) as handle:
            self._mutations = json.load(handle)

        with open(self._tool_inputs['JSON_pheno'][0].path) as handle:
            for row in json.load(handle):
                if row['category'] not in self._ab_info_by_category:
                    self._ab_info_by_category[row['category']] = []
                self._ab_info_by_category[row['category']].append(row)

    @staticmethod
    def __format_locus_name(str_in: str) -> str:
        """
        Formats the locus name.
        :param str_in: Input string
        :return: Formatted locus name
        """
        if str_in.endswith('_prom'):
            return str_in.replace('_prom', '<sub>prom</sub>')
        return str_in

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__parse_input_files()

        # Add content
        self.__add_resistance_type_overview()
        self.__add_antibiotics_table()

        # Add tables with mutations
        self.__add_mutations_table(
            'filtered', [m for m in self._mutations if m['passes_filt'] is True and m['lofreq'] is False],
            header='Mutations in resistance regions (filtered)')
        self.__add_mutations_table(
            'unfiltered', [m for m in self._mutations if m['passes_filt'] is False and m['lofreq'] is False],
            header='Mutations in resistance regions (unfiltered)')

        # Add table with low frequency mutations
        self.__add_mutations_table(
            'lofreq', [m for m in self._mutations if m['lofreq'] is True],
            header='Additional mutations in resistance regions detected by LoFreq')

        # Add visualization and DB info
        self.__add_visualization(self._tool_inputs['PNG'][0].path)
        self.__add_database_information()

        # Set tool output
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

        # Save JSON file in output directory
        self._section.add_file(self._tool_inputs['JSON'][0].path, Path(self._sub_folder, 'mutations.json'))

    def __add_antibiotics_table(self) -> None:
        """
        Adds an overview table with all antibiotics.
        :return: None
        """
        div = HtmlElement('div', attributes=[('class', 'border_bottom')])
        div.add_header('Predicted AMR', 4)

        confidence_groups = [
            {'name': 'Associated with R', 'confidence_levels': (
                amrutils.ConfidenceLevel.ASSOC_R, amrutils.ConfidenceLevel.ASSOC_R_int)},
            {'name': 'Not associated with R', 'confidence_levels': (
                amrutils.ConfidenceLevel.ASSOC_S, amrutils.ConfidenceLevel.ASSOC_S_int)},
            {'name': 'Uncertain significance', 'confidence_levels': (
                amrutils.ConfidenceLevel.UNKNOWN, amrutils.ConfidenceLevel.NOT_IN_DB)}
        ]

        # Add table heading
        table_data = [
            [HtmlElement('th', 'Antibiotic', [('rowspan', 2)]),
             HtmlElement('th', 'Abbreviation', [('rowspan', 2)]),
             HtmlElement('th', 'Predicted phenotype', [('rowspan', 2)]),
             HtmlElement('th', 'Mutations', [('colspan', 5)])],
            [HtmlElement('th', confidence_group['name']) for confidence_group in confidence_groups],
        ]

        for category, all_ab_data in self._ab_info_by_category.items():
            table_data.append([HtmlElement('th', category, [('colspan', 7)])])
            for ab_data in all_ab_data:
                color = 'red' if ab_data['phenotype'].startswith('R') else 'green'
                row = [ab_data['name'], ab_data['abbreviation'], HtmlTableCell(ab_data['phenotype'], color=color)]
                for confidence_group in confidence_groups:
                    mutations = list(itertools.chain(
                        *[ab_data['mutations'].get(level.value, []) for level in
                          confidence_group['confidence_levels']]))
                    row.append(', '.join([
                        f"{AMRReporter.__format_locus_name(m['region']['locus'])} ({m['name']})" for m in mutations])
                               if len(mutations) > 0 else '-')
                table_data.append(row)
        div.add_table(table_data, None, [('class', 'data')])

        # Save in tabular format
        relative_path = self.__save_table(table_data, None, 'amr-overview.tsv')
        div.add_link_to_file('Download (TSV)', relative_path)

        # Explanation
        div.add_paragraph(
            "This table contains the mutations related to the resistance to the antibiotics. Mutations that did not "
            "pass the variant filtering are indicated with a '*'. ")
        self._section.add_html_object(div)

    def __add_mutation_notes(self, div: HtmlElement, mutations: list[dict]) -> None:
        """
        Adds an overview with the notes for the mutations (if there are any).
        :param div: Section with mutation table
        :param mutations: List of mutations
        :return: None
        """
        # Group comments by index
        comment_by_idx_ab = {}
        for idx, m in enumerate(mutations):
            for association in m['associations']:
                if association['comment'] is None:
                    continue
                key = (idx, association['antibiotic_short'])
                comment_by_idx_ab[key] = association['comment']

        # Don't add anything if there are no notes
        if len(comment_by_idx_ab) == 0:
            return

        # Add table with notes
        div.add_header('Notes', 4)
        table_data = []
        for (idx, ab), comment in comment_by_idx_ab.items():
            table_data.append([mutations[idx]['name_full'], ab, comment])
        div.add_table(table_data, ['Mutation', 'Antibiotic', 'Note'], [('class', 'data')])

    def __add_mutations_table(self, suffix: str, mutations: list[dict[str, Any]], header: str) -> None:
        """
        Adds a mutation table.
        :param mutations: Mutation data
        :param suffix: Suffix for the table name (e.g. 'filtered')
        :param header: Header for the table
        :return: None
        """
        div = HtmlElement('div', attributes=[('class', 'border_bottom')])
        div.add_header(header, 4)

        table_data = []
        for row in mutations:

            # Skip synonymous mutations
            if row['effect'] == 'synonymous':
                logger.info(f"Skipping synonymous mutation: {row}")
                continue

            # Combine associations
            ab_str, conf_str = amrutils.combine_associations(row['associations'])

            table_data.append([
                AMRReporter.__format_locus_name(row['region']['locus']),
                row['position'],
                amrutils.shorten_nucleotide_str(row['ref']),
                amrutils.shorten_nucleotide_str(row['alt']),
                row['name'],
                ', '.join([str(x) for x in self._actg_counts_by_pos.get(row['position'], ['-'])]) if
                row['variant_type'] == 'snp' else 'n/a',
                ab_str,
                conf_str
            ])
        header = ['Locus', 'Position', 'Ref', 'Alt', 'Mutation', 'ACTG counts', 'Antibiotic', 'Level']

        if len(table_data) == 0:
            div.add_paragraph('No mutations found.')
            self._section.add_html_object(div)
            return

        # Add table
        div.add_table(table_data, header, [('class', 'data')])
        relative_path = self.__save_table(table_data, header, f'amr-variants-{suffix}.tsv')
        div.add_link_to_file('Download (TSV)', relative_path)
        self.__add_mutation_notes(div, mutations)
        self._section.add_html_object(div)

    def __add_resistance_type_overview(self) -> None:
        """
        Adds an overview with the resistance type for the sample.
        :return: None
        """
        with open(self._tool_inputs['JSON_amr_type'][0].path) as handle:
            data_amr_type = json.load(handle)
        div = HtmlElement('div', attributes=[('class', 'border_bottom')])
        div.add_header('Predicted drug susceptibility', 4)
        table_data = []
        for key, description in amrtypedetermination.DESCRIPTION_BY_RES_TYPE.items():
            table_data.append(
                ['☑' if data_amr_type['resistance_type'] == key else '☐', description]
            )
        div.add_table(table_data, table_attributes=[('class', 'information')])
        self._section.add_html_object(div)

    def __add_visualization(self, image_path: Path) -> None:
        """
        Adds the visualization of the mutations.
        :param image_path: Image path
        :return: None
        """
        div = HtmlElement('div', attributes=[('class', 'border_bottom')])
        div.add_header('Visualization', 4)
        relative_path = Path(self._sub_folder, 'resistance.png')
        self._section.add_file(image_path, relative_path)
        img = HtmlElement('img', attributes=[
            ('src', str(relative_path)), ('alt', 'visualization'), ('height', '960'), ('width', '960')])
        div.add_html_object(img)
        table_data = [
            [HtmlElement('th', 'Legend', [('colspan', '2')])],
            [HtmlTableCell('', color='plot_cds'), 'CDS'],
            [HtmlTableCell('', color='plot_rrna'), 'rRNA'],
            [HtmlTableCell('', color='plot_intergenic'), 'Intergenic'],
            [HtmlTableCell('', color='plot_pseudo'), 'pseudogene']
        ]
        div.add_table(table_data, table_attributes=[('class', 'data')])
        self._section.add_html_object(div)

    def __add_database_information(self) -> None:
        """
        Adds the database information to the report.
        :return: None
        """
        self._section.add_header('Database information', 4)
        self._section.add_paragraph(
            "Mutations are queried against the WHO catalogue of mutations associated with drug resistance in "
            f"<i>Mycobacterium tuberculosis</i> (version <b>{self._input_informs['screen']['version']}</b>), "
            "supplemented by a selection of mutations provided by the NRC. The following BED file contains the regions "
            "that are screened for resistance.")
        relative_path = self._section.add_file(self._tool_inputs['BED'][0].path, Path(
            self._sub_folder, 'regions_resistance.bed'))
        self._section.add_link_to_file('Resistance regions (BED)', relative_path)

    def __save_table(self, table_data: list[list[Any]], header: Optional[list[str]], filename: str) -> Path:
        """
        Saves the table to the report output.
        :param table_data: Table data
        :param header: Header
        :param filename: Filename
        :return: None
        """
        # Reformat columns with HTML tags
        for row in table_data:

            # Other HTML tags
            for i, cell in enumerate(row):
                if isinstance(cell, str):
                    in_str = cell
                elif isinstance(cell, HtmlElement):
                    in_str = cell.text
                else:
                    continue
                if '<sub>prom</sub>' in in_str:
                    in_str = in_str.replace('<sub>prom</sub>', '_prom')
                row[i] = re.sub('<[a-z]+>(.*?)</[a-z]+>', r'[\g<1>]', in_str)

        # Save in tabular format
        output_path = self.folder / filename
        TsvExporter.export(table_data, header, str(output_path))
        relative_path = Path(self._sub_folder) / output_path.name
        self._section.add_file(output_path, relative_path)
        return relative_path
