import json
from pathlib import Path
from typing import Any, Callable

from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlelement import HtmlElement
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltablecell import HtmlTableCell

from camel.app.core.tool import Tool


class SISTRReporter(Tool):
    """
    Creates an HTML report for the SISTR JSON output.
    """

    TITLE = 'SISTR'
    ANTIGEN_CONFIG = {
        'h1': {
            'loci': ['fliC'],
            'short': 'h1',
            'hits': lambda d: [d['h1_flic_prediction']],
        },
        'h2': {
            'loci': ['fljB'],
            'short': 'h2',
            'hits': lambda d: [d['h2_fljb_prediction']],
        },
        'o': {
            'loci': ['wzx', 'wzy'],
            'short': 'serogroup',
            'hits': lambda d: [
                d['serogroup_prediction']['wzx_prediction'],
                d['serogroup_prediction']['wzy_prediction'],
            ],
        },
    }
    CGMLST_TABLE = [
        {'key': 'cgmlst_ST', 'name': 'ST'},
        {'key': 'cgmlst_found_loci', 'name': 'Detected loci'},
        {'key': 'cgmlst_matching_alleles', 'name': 'Matching alleles'},
        {'key': 'cgmlst_subspecies', 'name': 'Subspecies', 'fmt': lambda x: f'<i>{x}</i>'},
        {'key': 'serovar_cgmlst', 'name': 'Serovar'},
    ]

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('SISTR Reporter', '0.1')
        self._section: HtmlReportSection | None = None

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Parse the input
        path_json_in = self._tool_inputs['JSON_SISTR'][0].path
        with path_json_in.open('r') as handle:
            data_sistr: dict[str, Any] = json.load(handle)[0]

        # Initialize the report section
        self._section = HtmlReportSection(
            SISTRReporter.TITLE,
            subtitle=self._input_informs['serotyping_sistr']['_name_full'],
        )

        # Add content
        self.__add_section_sistr(data_sistr)
        self._add_link_to_output(path_json_in)
        self._section.add_header('Additional information', 3)
        self._add_db_info(self._tool_inputs['DIR_sistr'][0].path)

        # Store the output
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_section_sistr(self, data_sistr: dict) -> None:
        """
        Adds the SISTR section to the HTML report.
        :param data_sistr: SISTR output data
        :return: None
        """
        if data_sistr['qc_status'] != 'PASS':
            self._section.add_warning_message(
                'SISTR QC failed: the isolate may be contaminated or have insufficient sequencing coverage.'
            )
        self._section.add_header('Output', 3)
        self._section.add_header('H1 (<i>fliC</i>)-type', 4)
        self._add_antigen_serotype_table_sistr(data_sistr, 'h1')
        self._section.add_header('H2 (<i>fljB</i>)-type', 4)
        self._add_antigen_serotype_table_sistr(data_sistr, 'h2')
        self._section.add_header('O-type', 4)
        self._add_antigen_serotype_table_sistr(data_sistr, 'o')
        self._add_cgmlst_results(data_sistr)
        self._section.add_header('Conclusion', 3)
        self._add_conclusion()

    @staticmethod
    def _fmt_value(val: Any, fmt: Callable | None = None) -> str:
        """
        Formats the given value.
        :param val: value to format as str
        :param: formatter function for formatting
        :return: Formatted value
        """
        if val is None:
            return 'n/a'
        if fmt is not None:
            return fmt(val)
        return str(val)

    def _add_cgmlst_results(self, data_sistr: dict) -> None:
        """
        Adds the cgMLST results.
        :param data_sistr: Parsed SISTR data
        :return: None
        """
        self._section.add_header('cgMLST', 3)
        self._section.add_table(
            [[SISTRReporter._fmt_value(data_sistr[r['key']], fmt=r.get('fmt')) for r in SISTRReporter.CGMLST_TABLE]],
            column_names=[r['name'] for r in SISTRReporter.CGMLST_TABLE],
            table_attributes=[('class', 'data')]
        )

    def _add_antigen_serotype_table_sistr(self, data_sistr: dict, antigen: str) -> None:
        """
        Generates and adds the table of an antigen type for sistr tool.
        :param data_sistr: SISTR output data
        :param: antigen, the salmonella antigen in lowercase, either o, h1, or h2
        :return: None
        """
        config = SISTRReporter.ANTIGEN_CONFIG.get(antigen)
        if not config:
            raise ValueError(f"Antigen must be one of: {', '.join(config.keys())}")

        best_hits_per_locus = config['hits'](data_sistr)
        locus_full = config['loci']
        locus_short = config['short']

        table_data = []
        header = [
            'Locus',
            '% Identity',
            'HSP/Locus length',
            'Contig',
            'Position in contig',
            'Predicted serotype',
        ]
        any_best_hits_found = False

        for i, best_hit in enumerate(best_hits_per_locus):
            if not best_hit.get('is_missing', True):
                any_best_hits_found = True
                top = best_hit.get('top_result', {})
                color = self._get_hit_color(best_hit)

                row = [
                    locus_full[i],
                    f"{top.get('pident', 0):.2f}",
                    f"{top.get('length', 0)}/{top.get('qlen', 0)}",
                    top.get('stitle', ''),
                    f"{top.get('sstart', '')}...{top.get('send', '')}",
                    best_hit.get(locus_short, ''),
                ]
                table_data.append([HtmlTableCell(x, color) for x in row])

            elif antigen == 'o':
                table_data.append(
                    [
                        HtmlTableCell(locus_full[i], 'red'),
                        HtmlElement(
                            'td',
                            'No match found',
                            attributes=[('colspan', 5), ('class', 'red')],
                        ),
                    ]
                )

        if any_best_hits_found:
            self._section.add_table(table_data, header, [('class', 'data')])
            if antigen == 'o':
                self._section.add_paragraph(
                    f"Predicted O antigen based on H antigens and serogroup: {data_sistr.get('o_antigen', 'N/A')}"
                )
        elif antigen in ('h1', 'h2'):
            self._section.add_paragraph(f"<i>{locus_full[0]}</i>: No match found")
        else:
            self._section.add_paragraph(
                f"<i>{locus_full[0]}</i> and <i>{locus_full[1]}</i>: No match found"
            )

    @staticmethod
    def _get_hit_color(hit: dict[str, Any]) -> str:
        """
        Returns the color for the corresponding hit.
        :param hit: Input hit
        :return: Color as a string
        """
        if (
            hit['top_result']['pident'] == 100.0
            and hit['top_result']['length'] == hit['top_result']['qlen']
        ):
            color = 'green'
        elif hit['top_result']['length'] == hit['top_result']['qlen']:
            color = 'lightgreen'
        else:
            color = 'grey'
        return color

    def _add_conclusion(self) -> None:
        """
        Generates and adds the conclusion sentence for sistr tool.
        :return: None
        """
        with self._tool_inputs['JSON_SISTR'][0].path.open('r') as handle:
            json_data = json.load(handle)[0]
        table_data = [
            (
                'Antigenic formula (O:H1:H2)',
                ':'.join(
                    [
                        str(json_data['o_antigen']),
                        str(json_data['h1']),
                        str(json_data['h2']),
                    ]
                ),
            ),
            ('Serotype', json_data['serovar']),
        ]
        self._section.add_table(
            table_data, ['Output', 'Prediction'], table_attributes=[('class', 'data')]
        )

    def _add_link_to_output(self, path_json: Path) -> None:
        """
        Adds a download link to the SISTR output file.
        :param path_json: JSON output file
        :return: None
        """
        relative_path = Path('serotyping', 'sistr', 'summary_out.json')
        self._section.add_link_to_file("Download (JSON)", relative_path)
        self._section.add_file(path_json, relative_path)

    def _add_db_info(self, db_dir: Path) -> None:
        """
        Adds the database information.
        :param db_dir: Input database directory
        :return: None
        """
        db_metadata_file = db_dir / 'db_update_info.json'
        if not db_metadata_file.is_file():
            raise FileNotFoundError(f'Database metadata not found: {db_metadata_file}')
        with db_metadata_file.open() as handle:
            metadata = json.load(handle)
            last_update_date = metadata['last_update_date']
        self._section.add_paragraph(f'Last database update: {last_update_date}')
