import json
from pathlib import Path
from typing import Any

from camel.app.core.reports.htmlelement import HtmlElement
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class SistrReporter(Tool):
    """
    Parses Sistr's JSON output results  and returns an HTML report.
    """
    TITLE = 'SISTR'

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Sistr Reporter', '0.1')
        self._section = None

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._section = HtmlReportSection(
            SistrReporter.TITLE, subtitle=self._input_informs['serotyping_sistr']['_name'])
        self.__add_section_sistr()
        self.__add_file_output()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_section_sistr(self) -> None:
        """
        Adds the SISTR section to the HTML report.
        :return: None
        """
        with self._tool_inputs['JSON_SISTR'][0].path.open('r') as handle:
            data_sistr = json.load(handle)[0]
            if data_sistr['qc_status'] == 'PASS':
                self._section.add_header('H1 (fliC)-type', 4)
                self.___add_antigen_serotype_table_sistr('h1')
                self._section.add_header('H2 (fljB)-type', 4)
                self.___add_antigen_serotype_table_sistr('h2')
                self._section.add_header('O-type', 4)
                self.___add_antigen_serotype_table_sistr('o')
                self._section.add_header('Conclusion', 4)
                self.___add_conclusion_sistr()
            else:
                self._section.add_paragraph('SISTR did not pass QC: the isolate might be contaminated, '
                                            'or not contain enough sequencing reads, or ...')
        db_dir = self._tool_inputs['DIR_sistr'][0].path
        self.___add_database_information(db_dir)

    def ___add_antigen_serotype_table_sistr(self, antigen: str) -> None:
        """
        Generates and adds the table of an antigen type for sistr tool.
        :param: antigen, the salmonella antigen in lowercase, either o, h1, or h2
        :return: None
        """
        table_data = []
        header = ['Locus', '% Identity', 'HSP/Locus length', 'Contig', 'Position in contig', 'Predicted serotype']
        with self._tool_inputs['JSON_SISTR'][0].path.open('r') as json_file:
            handle = json.load(json_file)[0]
        if antigen == 'h1':
            best_hits_per_locus = [handle['h1_flic_prediction']]
            locus_full = ['fliC']
            locus_short = ['h1']
        elif antigen == 'h2':
            best_hits_per_locus = [handle['h2_fljb_prediction']]
            locus_full = ['fljB']
            locus_short = ['h2']
        elif antigen == 'o':
            best_hits_per_locus = [handle['serogroup_prediction']['wzx_prediction'],
                                   handle['serogroup_prediction']['wzy_prediction']]
            locus_full = ['wzx', 'wzy']
            locus_short = ['serogroup']
        else:
            raise ValueError('antigen needs to be either h1, h2, or o !')

        any_best_hits_found = False
        # loop over best hits per locus; h1 and h2 will both only have one best hit each,
        # o has two loci, and therefore two best hits
        for index, best_hit in enumerate(best_hits_per_locus):
            if str(best_hit['is_missing']) == 'False':
                any_best_hits_found = True
                color = self._get_hit_color(best_hit)
                row = [
                    locus_full[index], format(best_hit['top_result']['pident'], '.2f'),
                    '/'.join([str(best_hit['top_result']['length']),
                              str(best_hit['top_result']['qlen'])]),
                    best_hit['top_result']['stitle'],
                    '...'.join([str(best_hit['top_result']['sstart']), str(best_hit['top_result']['send'])]),
                    best_hit[locus_short[0]]
                ]
                row = [HtmlTableCell(x, color) for x in row]
                table_data.append(row)
            else:
                if antigen == 'o':
                    row = [HtmlTableCell(locus_full[index], 'red'),
                           HtmlElement('td', 'No match found', attributes=[('colspan', 5), ('class', 'red')])]
                    table_data.append(row)

        if any_best_hits_found:
            self._section.add_table(table_data, header, [('class', 'data')])
            if antigen == 'o':
                self._section.add_paragraph(' '.join(['Predicted O antigen based on H antigens and serogroup:',
                                                      handle['o_antigen']]))
        else:  # if not any_best_hits_found:
            if antigen == 'h1' or antigen == 'h2':
                self._section.add_paragraph(f'{locus_full[0]}: No match found')
            else:
                self._section.add_paragraph(f'{locus_full[0]} and {locus_full[1]}: No match found')

    @staticmethod
    def _get_hit_color(hit: dict[str, Any]) -> str:
        """
        Gets a hits color based on its percentage identity and coverage.
        :return: color as str
        """
        if hit['top_result']['pident'] == 100.0 and \
                hit['top_result']['length'] == hit['top_result']['qlen']:
            color = 'green'
        elif hit['top_result']['length'] == hit['top_result']['qlen']:
            color = 'lightgreen'
        else:
            color = 'grey'
        return color

    def ___add_conclusion_sistr(self) -> None:
        """
        Generates and adds the conclusion sentence for sistr tool.
        :return: None
        """
        with self._tool_inputs['JSON_SISTR'][0].path.open('r') as handle:
            json_data = json.load(handle)[0]
        self._section.add_paragraph(' '.join(['Predicted antigenic profile (O:H1:H2):',
                                              ':'.join([str(json_data['o_antigen']),
                                                        str(json_data['h1']),
                                                        str(json_data['h2'])])]))
        self._section.add_paragraph(' '.join(['Predicted serotype:', json_data['serovar']]))

    def __add_file_output(self) -> None:
        """
        Add the output tsv file to the html.
        :return: None
        """
        relative_path = Path('serotyping', 'sistr', 'summary_out.json')
        self._section.add_link_to_file("Download (JSON)", relative_path)
        self._section.add_file(self._tool_inputs['JSON_SISTR'][0].path, relative_path)

    def ___add_database_information(self, db_dir: Path) -> None:
        """
        Adds the date of latest database update.
        :param db_dir: Input database directory
        :return: None
        """
        db_metadata_file = db_dir / 'db_update_info.json'
        if not db_metadata_file.is_file():
            raise FileNotFoundError(f'Database metadata not found: {db_metadata_file}')
        with db_metadata_file.open() as handle:
            metadata = json.load(handle)
            last_update_date = metadata['last_update_date']
        self._section.add_paragraph(f'Last updated: {last_update_date}')
