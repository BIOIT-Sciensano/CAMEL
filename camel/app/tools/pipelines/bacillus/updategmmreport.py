import ast

import pandas as pd
from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltablecell import HtmlTableCell

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class UpdateGMMReport(Tool):
    """
    Updates the GMM gene detection report with warning about detected GMMs.
    """
    INPUT_KEYS = ['TSV_STRAINS', 'TSV_GMM_VECTORS', 'TSV_GMM_JUNCTIONS', 'VAL_HTML_VECTORS', 'VAL_HTML_JUNCTIONS',
                  'TSV_GMM_DB']
    COLOR_CODE = {'STRAIN_MATCH': 'green', 'GMM_MATCH': 'yellow', 'BOTH_MATCH': 'red', 'UNKNOWN_MATCH': 'grey'}

    def __init__(self) -> None:
        """
        Initializes the tool.
        :return: None
        """
        super().__init__('UpdateGMMReport', '0.1')

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if any(key not in self._tool_inputs for key in self.INPUT_KEYS):
            raise InvalidToolInputError(
                "Tool requires {} inputs".format(', '.join(UpdateGMMReport.INPUT_KEYS)))
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the tool.
        :return: None
        """
        output_report = self._update_report()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(output_report), self._tool_inputs['VAL_HTML_JUNCTIONS'][0]]

    def _update_report(self) -> HtmlReportSection:
        """
        Updates the report with the GMM warning table.
        :return: Updated report section
        """
        matches = self._parse_tsv_files()
        colored_matches = self._generate_table_with_colors(matches)
        current_report_section = self._tool_inputs['VAL_HTML_VECTORS'][0].value

        new_section = HtmlReportSection('Interpretation')

        if not (colored_matches['raw_table']['strain'] and colored_matches['raw_table']['construct']):
            new_section.add_paragraph('No GMM construct detected.')
        else:
            column_names = ['strain', 'construct']
            new_section.add_table(colored_matches['colored_table'], column_names, [('class', 'data')])

            if matches['strain']:
                new_section.add_paragraph(f'The strain matches closely to '
                                          f'<b>{colored_matches["raw_table"]["strain"]}</b> '
                                          f'which has been used in GMMs.')
            else:
                new_section.add_paragraph('The strain does not match any known GMM strains in the database.')
            if matches['construct']:
                new_section.add_paragraph(f'The <b>{colored_matches["raw_table"]["construct"]}</b> '
                                          f'transgenic construct was detected in the strain.')
            else:
                new_section.add_paragraph(
                    'No transgenic constructs from the database were detected in the strain.')
            self.__add_explanation_matches(new_section)

        new_section.add_warning_message('The pipeline uses a targeted approach, which means that constructs '
                                        'and/or strains that are not in the database will be missed.')
        current_report_section.add_horizontal_line()
        current_report_section.add_html_object(new_section)
        return current_report_section

    def _parse_tsv_files(self) -> dict:
        """
        Parses the TSV files passed as input.
        :return: Dictionary with match, or False if no match is found
        """
        strain_hits = []
        for f in self._tool_inputs['TSV_STRAINS']:
            with open(f.path) as handle:
                tsv_strain = pd.read_csv(handle, sep='\t', names=['metric', 'value'])
                rows_of_interest = [tsv_strain['metric'].tolist().index(x) for x in tsv_strain['metric'] if
                                    'closest_strain' in x]
                if len(rows_of_interest) == 1:
                    strain_hits.append(tsv_strain.loc[rows_of_interest[0], 'value'])
                else:
                    strain_hits.extend(tsv_strain.loc[rows_of_interest, 'value'].tolist())
        gmm_hits = []
        with open(self._tool_inputs['TSV_GMM_VECTORS'][0].path) as handle:
            all_lines = handle.readlines()
            gmm_hits_list = ast.literal_eval(all_lines[0].strip().split('\t')[1])
            gmm_hits.extend(hit[1] for hit in gmm_hits_list)

        return {'strain': strain_hits,
                'construct': gmm_hits}

    def _generate_table_with_colors(self, match_table: dict) -> dict:
        """
        Generates a table with rows colored.
        :param match_table: List of table to be colored
        :return: Colored table
        """
        tsv_gmm_db = pd.read_csv(self._tool_inputs['TSV_GMM_DB'][0].path).to_dict(orient='records')
        perfect_matches = [(hit['construct'], hit['strain']) for hit in tsv_gmm_db]
        identified_constructs = [x['construct'] for x in tsv_gmm_db]
        identified_strains = [x['strain'] for x in tsv_gmm_db]

        strain_hit = [s for s in match_table['strain'] if s in identified_strains]
        gmm_hit = [s for s in match_table['construct'] if s in identified_constructs]

        if not (strain_hit and gmm_hit):
            return {'raw_table': {'strain': [], 'construct': gmm_hit}}

        color = None
        if strain_hit and not gmm_hit:
            color = UpdateGMMReport.COLOR_CODE['STRAIN_MATCH']
            gmm_hit = [None]
        if gmm_hit and not strain_hit:
            color = UpdateGMMReport.COLOR_CODE['GMM_MATCH']
            strain_hit = [None]
        if strain_hit and gmm_hit and (*gmm_hit, *strain_hit) in perfect_matches:
            color = UpdateGMMReport.COLOR_CODE['BOTH_MATCH']
        if strain_hit and gmm_hit and (*gmm_hit, *strain_hit) not in perfect_matches:
            color = UpdateGMMReport.COLOR_CODE['UNKNOWN_MATCH']
        temp_table = [HtmlTableCell(text, color=color) for text in [*strain_hit, *gmm_hit]]
        table_to_return = {'colored_table': [(temp_table[0], temp_table[1])],
                           'raw_table': {'strain': strain_hit[0], 'construct': gmm_hit[0]}}
        return table_to_return

    def __add_explanation_matches(self, section: HtmlReportSection) -> None:
        """
        Adds information about the different type of matches to the bottom of the report.
        :param section: Report section
        :return: None
        """
        section.add_header('Extra information', 3)
        section.add_paragraph('The following colors are used to denote the different type of hits:')
        section.add_table([
            [HtmlTableCell('', color=UpdateGMMReport.COLOR_CODE['STRAIN_MATCH']), 'Matching frequent GMM strain'],
            [HtmlTableCell('', color=UpdateGMMReport.COLOR_CODE['GMM_MATCH']), 'Matching GMM construct'],
            [HtmlTableCell('', color=UpdateGMMReport.COLOR_CODE['BOTH_MATCH']),
             'Both GMM strain and construct detected'],
            [HtmlTableCell('', color=UpdateGMMReport.COLOR_CODE['UNKNOWN_MATCH']),
             'GMM strain and construct detected - Unknown combination'],
        ], None, [('class', 'data')])
