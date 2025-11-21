import re
from pathlib import Path

import pandas as pd

from camel.app.core.reports.htmlelement import HtmlElement
from camel.app.core.reports.htmlexpandablediv import HtmlExpandableDiv
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class MOBReconReporter(Tool):
    """
    Creates an HTML report for MOB-recon.
    """

    COLUMN_MAPPING = {
        'id': {'title': 'ID'},
        'num_contigs': {'title': 'Nb. of contigs'},
        'size': {'title': 'Size', 'fmt': lambda x: f'{x:,}'},
        'gc': {'title': '% GC-content', 'fmt': lambda x: f'{x * 100:.2f}'},
        'predicted_mobility': {'title': 'Pred. mobility'},
        'rep_type(s)': {'title': 'Rep. types', 'fmt': lambda x: x.replace(',', ', ')},
        'relaxase_type(s)': {'title': 'Relaxase types', 'fmt': lambda x: x.replace(',', ', ')}
    }

    CONTIG_COLUMN_MAPPING = {
        'contig_id': {'title': 'ID'},
        'molecule_type': {'title': 'Mol. type'},
        'primary_cluster_id': {'title': 'Prim. cluster ID'},
        'secondary_cluster_id': {'title': 'Sec. cluster ID'},
        'size': {'title': 'Size', 'fmt': lambda x: f'{x:,}'},
        'gc': {'title': '% GC-content', 'fmt': lambda x: f'{x * 100:.2f}'}
    }

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('MOB-recon reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError('TSV input is required')
        if 'TSV_contigs' not in self._tool_inputs:
            raise InvalidToolInputError('TSV contigs input is required')
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA input is required')
        if 'mob_recon' not in self._input_informs:
            raise InvalidToolInputError("MOB-recon informs input is required")
        super()._check_input()

    def _add_overview_table(self, section: HtmlReportSection) -> None:
        """
        Adds an overview table to the report.
        :param section: Report section
        :return: None
        """
        if len(self._input_informs['mob_recon']['detected_plasmids']) > 0:
            data_overview = pd.read_table(self._tool_inputs['TSV'][0].path)
        else:
            section.add_paragraph('No plasmids detected')
            return

        data_overview['id'] = data_overview['sample_id'].apply(lambda x: re.search('.*:(.*)', x).group(1))
        data_overview['id'] = data_overview['id'].apply(lambda x: MOBReconReporter.format_plasmid_id(x))
        section.add_header('Overview', 3)
        table_data = [
            [d.get('fmt', lambda x: x)(row[col]) for
             col, d in MOBReconReporter.COLUMN_MAPPING.items()] for row in data_overview.to_dict('records')
        ]

        # Add column with a download link for FASTA files
        for row in table_data:
            id_ = row[0]
            path_fasta = next(io.path for io in self._tool_inputs['FASTA'] if id_ in io.path.name)
            relative_path = Path('mob-suite', path_fasta.name)
            section.add_file(path_fasta, relative_path)
            # noinspection PyTypeChecker
            row.append(HtmlTableCell('Download (FASTA)', link=str(relative_path)))

        # Add table
        section.add_table(
            table_data, [c['title'] for c in MOBReconReporter.COLUMN_MAPPING.values()] + ['Sequence'],
            [('class', 'data')])

    def _add_contig_report(self, section: HtmlReportSection) -> None:
        """
        Adds a contig report to the report.
        :param section: Report section
        :return: None
        """
        contig_data = pd.read_table(self._tool_inputs['TSV_contigs'][0].path,
                                    usecols=list(MOBReconReporter.CONTIG_COLUMN_MAPPING.keys()))
        reordered_contig_data = contig_data[list(MOBReconReporter.CONTIG_COLUMN_MAPPING.keys())]
        table_data = [
            [d.get('fmt', lambda x: x)(row[col])
             for col, d in MOBReconReporter.CONTIG_COLUMN_MAPPING.items()]
            for row in reordered_contig_data.to_dict('records')
        ]

        if len(table_data) > 10:
            div = HtmlExpandableDiv("contig-overview", f'{len(table_data)} rows.')
        else:
            div = HtmlElement('div')

        section.add_header('Contig overview', 3)
        div.add_table(table_data, [c['title'] for c in MOBReconReporter.CONTIG_COLUMN_MAPPING.values()],
                      [('class', 'data')])
        section.add_html_object(div)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('MOB-recon', subtitle=self._input_informs['mob_recon']['_name_full'])
        self._add_overview_table(section)
        if 'contig_report' in self._parameters:
            self._add_contig_report(section)

        # Download overview
        relative_path = Path('mob-suite', 'mob_recon.tsv')
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path)
        section.add_link_to_file('Download complete output (TSV)', relative_path)

        # Download contig report
        relative_path = Path('mob-suite', 'contig_report.tsv')
        section.add_file(self._tool_inputs['TSV_contigs'][0].path, relative_path)
        section.add_link_to_file('Download contig report (TSV)', relative_path)

        # Tool output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]

    @staticmethod
    def format_plasmid_id(str_in: str) -> str:
        """
        Formats the plasmid id (reduces the length for novel plasmid clusters).
        :param str_in: Input string
        :return: Formatted string
        """
        m = re.match(r'novel_(\w+)', str_in)
        if not m:
            return str_in
        return f'novel_{m.group(1)[:4]}'
