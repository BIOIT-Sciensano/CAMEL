from pathlib import Path

import pandas as pd
from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlexpandablediv import HtmlExpandableDiv
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltablecell import HtmlTableCell
from camelcore.app.reports.htmltableformatter import FormatEntry, HtmlTableFormatter

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class ReporterDbComparison(Tool):
    """
    Creates an output report for the database comparison.
    """

    COLS_MATCHES: list[FormatEntry] = [
        {'key': 'rank', 'title': 'Rank'},
        {'key': 'accession', 'title': 'Accession'},
        {'key': 'identity', 'title': '% Identity', 'fmt': HtmlTableFormatter.FLOAT_FMT},
        {'key': 'differences', 'title': 'Differences', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'mismatches', 'title': 'Mismatches', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'indels', 'title': 'Indels', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'coverage', 'title': '% Aligned', 'fmt': lambda x: f'{100 * x:.1f}'},
    ]
    COLS_METADATA = ['Host', 'Country', 'Collection_Date']

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Reporter: Database comparison', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        toolutils.check_input(self, keys_required=['TSV', 'DB'])
        if 'db_comparison' not in self._input_informs:
            raise InvalidToolInputError('Database comparison informs are required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        """
        data_matches = pd.read_table(self._tool_inputs['TSV'][0].path, keep_default_na=False, na_values=['NA', ''])
        section = HtmlReportSection('Database comparison')
        section.add_table([
            ['Database:', self._tool_inputs['DB'][0].path.name]
        ], table_attributes=[('class', 'information')])
        self.__add_overview_table(section, data_matches)
        self.__add_matches_tables(section, data_matches)
        self.__add_output_files(section)
        self.__add_information(section)
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]

    def __create_identity_cell(self, identity: float) -> HtmlTableCell:
        """
        Returns a table cell with the identity, colored based on the thresholds.
        :param identity: Identity (%)
        :return: Table cell
        """
        if identity >= float(self.get_param_value('identity_warn')):
            color = 'green'
        elif identity >= float(self.get_param_value('identity_fail')):
            color = 'yellow'
        else:
            color = 'red'
        return HtmlTableCell(f'{identity:.2f}', color=color)

    def __add_overview_table(self, section: HtmlReportSection, data_matches: pd.DataFrame) -> None:
        """
        Adds a table with the closest database match for each segment.
        :param section: Report section
        :param data_matches: Database comparison matches
        :return: None
        """
        section.add_header('Overview', 3)
        cols_meta = [c for c in ReporterDbComparison.COLS_METADATA if c in data_matches.columns]
        table_data = []
        for segment, data_seg in data_matches.groupby('segment'):
            closest = data_seg.iloc[0]
            table_data.append([
                segment,
                closest['accession'],
                self.__create_identity_cell(closest['identity']),
                f"{closest['differences']:,}",
                f"{100 * closest['coverage']:.1f}",
                f"{(data_seg['rank'] == 1).sum():,}",
                *[closest[c] if pd.notna(closest[c]) else 'n/a' for c in cols_meta]
            ])
        header = [
            'Segment', 'Closest match', '% Identity', 'Differences', '% Aligned', 'Nb. equally close matches',
            *[c.replace('_', ' ').capitalize() for c in cols_meta]]
        section.add_table(table_data, header, [('class', 'data')])

        # Warn if no sufficiently similar database match was found
        segments_low = [
            s for s, d in data_matches.groupby('segment')
            if d['identity'].iloc[0] < float(self.get_param_value('identity_fail'))
        ]
        if len(segments_low) > 0:
            segments_low_str = ','.join(str(x) for x in segments_low)
            section.add_alert(
                f"No closely related database sequence was found for segment(s) {segments_low_str} "
                f"(identity < {float(self.get_param_value('identity_fail')):.0f}%). "
                "Interpret the associated metadata with caution.",
                'warning')

    def __add_matches_tables(self, section: HtmlReportSection, data_matches: pd.DataFrame) -> None:
        """
        Adds a table with the top matches for each segment.
        :param section: Report section
        :param data_matches: Database comparison matches
        :return: None
        """
        section.add_header('Top database matches by segment', 3)
        div_ = HtmlExpandableDiv('db_comparison_matches', 'Matches by segment')
        for segment, data_seg in data_matches.groupby('segment'):
            div_.add_header(f'Segment {segment}', 3)
            div_.add_table(
                HtmlTableFormatter.format_table_data(data_seg, ReporterDbComparison.COLS_MATCHES),
                [c['title'] for c in ReporterDbComparison.COLS_MATCHES],
                [('class', 'data')])
        section.add_html_object(div_)

    def __add_output_files(self, section: HtmlReportSection) -> None:
        """
        Adds the output files.
        :param section: Report section
        :return: None
        """
        section.add_header('Output files', 3)
        path_rel = Path('db_comparison', 'closest_matches.tsv')
        section.add_file(self._tool_inputs['TSV'][0].path, path_rel)
        section.add_table([
            ['Closest matches (with metadata)', HtmlTableCell('Download (TSV)', link=str(path_rel))]
        ], ['File', 'Download'], [('class', 'data')])

    def __add_information(self, section: HtmlReportSection) -> None:
        """
        Adds the additional information.
        :param section: Report section
        :return: None
        """
        section.add_header('Additional information', 3)
        section.add_paragraph(
            'The final consensus sequence of each segment is compared to all sequences in the database. '
            'Candidate matches are identified using <i>mash</i> and aligned to the consensus sequence using '
            '<i>minimap2</i>.'
            'Matches are ranked by the number of differences (mismatches + indel events). '
            'Ambiguous positions are not compared and overhanging ends are not penalized. '
            'Matches covering less than <b>{:.0f}%</b> of the consensus are not reported. '.format(
                100 * float(self._input_informs['db_comparison']['min_coverage']),
            ))
