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
from camel.app.toolkits.viralconsensus import dbutils


class ReporterDbComparison(Tool):
    """
    Creates an output report for the database comparison.
    """

    COLS_MATCHES: list[FormatEntry] = [
        {'key': 'rank', 'title': 'Rank'},
        {'key': 'accession', 'title': 'Accession'},
        {'key': 'identity', 'title': 'Identity', 'fmt': lambda x: f'{x:.2f}%'},
        {'key': 'coverage', 'title': 'Aligned', 'fmt': lambda x: f'{100 * x:.2f}%'},
        {'key': 'differences', 'title': 'Differences', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'mismatches', 'title': 'Mismatches', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'indels', 'title': 'Indels', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'nb_accessions', 'title': 'Nb. accessions', 'fmt': lambda x: f'{int(x):,}' if pd.notna(x) else 'n/a'},
    ]
    COLS_METADATA = ['Genotype', 'Subtype', 'Host', 'Country', 'Collection_Date']
    TITLES_SUMMARY = {
        'Genotype': 'Genotypes', 'Subtype': 'Subtypes', 'Host': 'Hosts', 'Country': 'Countries',
        'Collection_Date': 'Collection dates'}
    NB_TOP_VALUES = 3

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Reporter: Database comparison', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        toolutils.check_input(self, keys_required=['TSV', 'TSV_equal', 'DB'])
        if 'db_comparison' not in self._input_informs:
            raise InvalidToolInputError('Database comparison informs are required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        """
        data_matches, data_equal = [pd.read_table(
            self._tool_inputs[key][0].path, keep_default_na=False, na_values=[''],
            dtype={c: str for c in ReporterDbComparison.COLS_METADATA}) for key in ('TSV', 'TSV_equal')]
        section = HtmlReportSection('Database comparison')
        manifest = dbutils.read_manifest(self._tool_inputs['DB'][0].path)
        section.add_table([
            ['Database:', manifest.get('database', self._tool_inputs['DB'][0].path.name)],
            ['Download date:', manifest.get('ncbi_download_date', 'n/a')]
        ], table_attributes=[('class', 'information')])
        if len(data_matches) > 0:
            self.__add_overview_table(section, data_matches, data_equal)
            self.__add_alerts(section, data_matches)
            self.__add_matches_tables(section, data_matches)
        else:
            section.add_alert('No matching database sequences found for the consensus sequence.', 'warning')
        self.__add_output_files(section)
        self.__add_information(section)
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]

    def __create_colored_cell(self, value: float, key: str, fmt: str) -> HtmlTableCell:
        """
        Returns a table cell with the value, colored based on the warning / fail thresholds.
        :param value: Value (%)
        :param key: Threshold parameter prefix ('identity' or 'coverage')
        :param fmt: Format string for the value
        :return: Table cell
        """
        if value >= float(self.get_param_value(f'{key}_warn')):
            color = 'green'
        elif value >= float(self.get_param_value(f'{key}_fail')):
            color = 'yellow'
        else:
            color = 'red'
        return HtmlTableCell(fmt.format(value), color=color)

    @staticmethod
    def _summarize_counts(values: pd.Series) -> str:
        """
        Summarizes the values as counts of the most common values, e.g. 'USA (10), New Zealand (4), other (3)'.
        :param values: Values
        :return: Summary
        """
        counts = values.fillna('n/a').value_counts().sort_index().sort_values(ascending=False, kind='stable')
        parts = [f'{value} ({count:,})' for value, count in counts.iloc[:ReporterDbComparison.NB_TOP_VALUES].items()]
        if len(counts) > ReporterDbComparison.NB_TOP_VALUES:
            parts.append(f'other ({counts.iloc[ReporterDbComparison.NB_TOP_VALUES:].sum():,})')
        return ', '.join(parts)

    @staticmethod
    def _summarize_dates(values: pd.Series) -> str:
        """
        Summarizes the dates as a range (year and month), e.g. '2018-10 - 2019-03'.
        :param values: Dates (YYYY, YYYY-MM or YYYY-MM-DD)
        :return: Summary
        """
        dates = values.dropna().sort_values()
        if len(dates) == 0:
            return 'n/a'
        date_first, date_last = dates.iloc[0][:7], dates.iloc[-1][:7]
        return date_first if date_first == date_last else f'{date_first} - {date_last}'

    def __add_overview_table(
            self, section: HtmlReportSection, data_matches: pd.DataFrame, data_equal: pd.DataFrame) -> None:
        """
        Adds a table with the closest database match for each segment, and a summary of the metadata of all equally
        close accessions.
        :param section: Report section
        :param data_matches: Database comparison matches
        :param data_equal: Equally close accessions
        :return: None
        """
        section.add_header('Overview', 3)
        cols_meta = [c for c in ReporterDbComparison.COLS_METADATA if c in data_equal.columns]
        table_data = []
        segments_db = dbutils.get_segments(self._tool_inputs['DB'][0].path)
        for segment in sorted(set(segments_db) | set(data_matches['segment'])):
            data_seg = data_matches[data_matches['segment'] == segment]

            # Segment without result
            if len(data_seg) == 0:
                table_data.append([
                    segment, 'n/a', HtmlTableCell('n/a', color='red'),
                    self.__create_colored_cell(0.0, 'coverage', '{:.2f}%'), 'n/a', 'n/a', 'n/a',
                    *['n/a' for _ in cols_meta]])
                continue

            # Closest match
            closest = data_seg.iloc[0]
            data_rank_1 = data_seg[data_seg['rank'] == 1]
            data_equal_seg = data_equal[data_equal['segment'] == segment]
            table_data.append([
                segment,
                closest['accession'],
                self.__create_colored_cell(closest['identity'], 'identity', '{:.2f}%'),
                self.__create_colored_cell(100 * closest['coverage'], 'coverage', '{:.2f}%'),
                f"{closest['differences']:,}",
                f"{len(data_rank_1):,}",
                f"{int(data_rank_1['nb_accessions'].sum()):,}" if data_rank_1['nb_accessions'].notna().all() else 'n/a',
                *[ReporterDbComparison._summarize_dates(data_equal_seg[c]) if c == 'Collection_Date' else
                  ReporterDbComparison._summarize_counts(data_equal_seg[c]) for c in cols_meta]
            ])
        header = [
            'Segment', 'Closest match', 'Identity', 'Aligned', 'Differences', 'Nb. equally close matches',
            'Nb. equally close accessions', *[ReporterDbComparison.TITLES_SUMMARY[c] for c in cols_meta]]
        section.add_table(table_data, header, [('class', 'data')])
        section.add_paragraph('The metadata is summarized for all equally close accessions.')

    def __add_alerts(self, section: HtmlReportSection, data_matches: pd.DataFrame) -> None:
        """
        Adds alerts for segments without result, closest matches with a low coverage or a low identity.
        :param section: Report section
        :param data_matches: Database comparison matches
        :return: None
        """
        # Warn for segments in the database without result
        segments_missing = [
            s for s in dbutils.get_segments(self._tool_inputs['DB'][0].path) if s not in set(data_matches['segment'])]
        if len(segments_missing) > 0:
            section.add_alert(
                f"No database comparison for segment(s) {', '.join(segments_missing)}: segment not detected, consensus "
                f"sequence without called bases or no matching database sequence.",
                'warning')

        # Warn if the closest match does not cover the consensus sequence
        coverage_fail = float(self.get_param_value('coverage_fail'))
        segments_low_cov = [
            s for s, d in data_matches.groupby('segment') if 100 * d['coverage'].iloc[0] < coverage_fail]
        if len(segments_low_cov) > 0:
            section.add_alert(
                f"The closest match covers less than {coverage_fail:.0f}% of the consensus sequence for segment(s) "
                f"{', '.join(segments_low_cov)} (incomplete consensus sequence with N positions, or shorter database "
                f"sequence). The identity and differences only apply to the aligned part, interpret the results with "
                f"caution.",
                'warning')

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
        cols_format = ReporterDbComparison.COLS_MATCHES + [
            {'key': c, 'title': c.replace('_', ' ').capitalize(), 'fmt': lambda x: x if pd.notna(x) else 'n/a'}
            for c in ReporterDbComparison.COLS_METADATA if c in data_matches.columns]
        div_ = HtmlExpandableDiv('db_comparison_matches', 'Matches by segment')
        for segment, data_seg in data_matches.groupby('segment'):
            div_.add_header(f'Segment {segment}', 3)
            div_.add_table(
                HtmlTableFormatter.format_table_data(data_seg, cols_format),
                [c['title'] for c in cols_format],
                [('class', 'data')])
        section.add_html_object(div_)

    def __add_output_files(self, section: HtmlReportSection) -> None:
        """
        Adds the output files.
        :param section: Report section
        :return: None
        """
        section.add_header('Output files', 3)
        table_data = []
        for key, description in (
                ('TSV', 'Closest matches (with metadata)'),
                ('TSV_equal', 'Equally close sequences, all accessions (with metadata)')):
            path_rel = Path('db_comparison', self._tool_inputs[key][0].path.name)
            section.add_file(self._tool_inputs[key][0].path, path_rel)
            table_data.append([description, HtmlTableCell('Download (TSV)', link=str(path_rel))])
        section.add_table(table_data, ['File', 'Download'], [('class', 'data')])

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
            '<i>minimap2</i>. '
            'Matches are ranked by the number of differences (mismatches + indel events). '
            'Ambiguous positions are not compared and overhanging ends are not penalized. '
            'Matches covering less than <b>{:.0f}%</b> of the called (non-N) positions of the consensus are not '
            'reported. The aligned fraction is relative to the full consensus sequence, including N positions. '.format(
                100 * float(self._input_informs['db_comparison']['min_coverage']),
            ))

        # Cutoffs used for the coloring
        div_ = HtmlExpandableDiv('db_comparison_color_coding', 'Color coding')
        identity_warn = float(self.get_param_value('identity_warn'))
        identity_fail = float(self.get_param_value('identity_fail'))
        div_.add_table([
            [HtmlTableCell(f'≥ {identity_warn:.0f}%', color='green'), 'Closely related database sequence'],
            [HtmlTableCell(f'{identity_fail:.0f} - {identity_warn:.0f}%', color='yellow'), 'Related database sequence'],
            [HtmlTableCell(f'< {identity_fail:.0f}%', color='red'),
             'No closely related database sequence, interpret the metadata with caution'],
        ], ['Identity (closest match)', 'Interpretation'], [('class', 'data')])

        coverage_warn = float(self.get_param_value('coverage_warn'))
        coverage_fail = float(self.get_param_value('coverage_fail'))
        div_.add_table([
            [HtmlTableCell(f'≥ {coverage_warn:.0f}%', color='green'), 'Consensus sequence covered by the match'],
            [HtmlTableCell(f'{coverage_fail:.0f} - {coverage_warn:.0f}%', color='yellow'),
             'Consensus sequence partially covered (e.g. shorter database sequence without UTRs)'],
            [HtmlTableCell(f'< {coverage_fail:.0f}%', color='red'),
             'Consensus sequence not covered (e.g. incomplete consensus sequence with N positions), interpret the '
             'results with caution'],
        ], ['Aligned (closest match)', 'Interpretation'], [('class', 'data')])
        section.add_html_object(div_)
