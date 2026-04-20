import json
from pathlib import Path
from typing import Any

from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.reports.htmlelement import HtmlElement
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.tool import Tool
from camel.app.loggers import logger


class SeqSero2Reporter(Tool):
    """
    Creates an output report for SeqSero2.
    """

    TITLE = 'SeqSero2'

    COLS = [
        {'name': 'O-antigen', 'key': 'O antigen prediction', 'name_short': 'o_antigen'},
        {'name': 'H1-antigen (<i>fliC</i>)', 'key': 'H1 antigen prediction(fliC)', 'name_short': 'h1_antigen'},
        {'name': 'H2-antigen (<i>fljB</i>)', 'key': 'H2 antigen prediction(fljB)', 'name_short': 'h2_antigen'},
        {'name': 'Antigenic formula', 'key': 'Predicted antigenic profile', 'name_short': 'formula'},
        {'name': 'Serotype', 'key': 'Predicted serotype', 'name_short': 'serotype'},
    ]

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('SeqSero2 Reporter', '0.1')
        self._section = None

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._section = HtmlReportSection(
            SeqSero2Reporter.TITLE,
            subtitle=self._input_informs['serotyping_seqsero2']['_name_full'],
        )
        self.__add_section_seqsero()
        self.__add_db_info()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_section_seqsero(self) -> None:
        """
        Adds the SeqSero2 section to the HTML report.
        :return: None
        """
        self._section.add_header('Output', 3)
        table_rows = []
        notes = []
        warnings = []

        # Potential output files from SeqSero2
        input_configs = [
            ('TXT_seqsero2_kmer', 'Assembly', 'Kmer'),
            ('TXT_seqsero2_allele', 'Reads', 'Allele'),
            ('TXT_seqsero2_kmerread', 'Reads', 'Kmer'),
        ]

        for key, label, mode in input_configs:
            informs_key = key.replace('TXT_', '')
            if key in self._tool_inputs:
                output = SeqSero2Reporter.parse_seqsero_output(self._tool_inputs[key][0].path)
                row_data = [output.get(k['key'], '') for k in self.COLS]
                table_rows.append([label, mode, *row_data])
                notes.append([label, mode, output.get('Note', 'n/a')])
                if self._input_informs.get(informs_key, {}).get('simulated'):
                    warnings.append(f"Input for mode '{mode}' on '{label}' was simulated from the assembly.")
            else:
                error_cell = HtmlElement(
                    'td', 'Not supported', [('colspan', len(self.COLS))]
                )
                table_rows.append([label, mode, error_cell])

        # Add main output table
        column_names = ['Input', 'Mode'] + [c['name'] for c in self.COLS]
        self._section.add_table(
            table_rows, column_names=column_names, table_attributes=[('class', 'data')]
        )

        # Add notes
        self._section.add_header('Notes', 4)
        self._section.add_table(
            notes, ['Input', 'Mode', 'Note'], table_attributes=[('class', 'data')]
        )

        # Add warnings (if applicable)
        for warning in warnings:
            self._section.add_warning_message(warning)

    @staticmethod
    def parse_seqsero_output(input_file_path: Path) -> dict[str, Any]:
        """
        Parses the SeqSero2 output file.
        :param input_file_path: Input file path
        :return: Parsed info as a dictionary
        """
        logger.debug(f'Parsing: {input_file_path}')
        output = {}
        with input_file_path.open('r') as handle:
            for line in handle:
                parts = line.rstrip().split(':\t')
                output[parts[0]] = parts[1].strip() if len(parts) > 1 else ''
        return output

    def __add_db_info(self) -> None:
        """
        Adds the database information.
        :return: None
        """
        self._section.add_header('Additional information', 3)
        db_dir = self._tool_inputs['DIR_seqsero2'][0].path
        db_metadata_file = db_dir / 'db_update_info.json'
        if not db_metadata_file.is_file():
            raise FileNotFoundError(f'Database metadata not found: {db_metadata_file}')
        with db_metadata_file.open() as handle:
            metadata = json.load(handle)
            last_update_date = metadata['last_update_date']
        self._section.add_paragraph(f'Last database update: {last_update_date}')
