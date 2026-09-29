import json
from pathlib import Path

import pandas as pd
from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlexpandablediv import HtmlExpandableDiv
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltablecell import HtmlTableCell
from camelcore.app.reports.htmltableformatter import FormatEntry, HtmlTableFormatter

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool
from camel.app.toolkits.viralconsensus import dbutils, seqidutils


class ReporterRefSelection(Tool):
    """
    Class to generate reports for the reference selection workflow.
    """

    COLS_SEGMENT: list[FormatEntry] = [
        {'key': 'ref_id_fmt', 'title': 'Reference'},
        {'key': 'length', 'title': 'Length', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'identity', 'title': '% Identity'}, # Already formatted
        {'key': 'p_val', 'title': 'P-value'},
        {'key': 'median_mult', 'title': 'Median-multiplicity (Cov.)', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'hashes', 'title': 'Matching hashes'},
        {'key': 'score_final', 'title': 'Score', 'fmt': HtmlTableFormatter.FLOAT_FMT},
    ]

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Reporter: Reference selection', '0.1')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Parse input
        with open(self._tool_inputs['JSON'][0].path) as handle:
            ref_info_by_seg = json.load(handle)
        data_mash = pd.read_table(self._tool_inputs['TSV'][0].path, keep_default_na=False, na_values='-')

        # Create the HTML report
        section = HtmlReportSection('Reference selection')
        self.__add_overview_table(section, ref_info_by_seg)
        self._add_table_ref_metadata(section, self._tool_inputs['DB'][0].path, ref_info_by_seg)
        self._add_segments_tables(section, data_mash)
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        if 'JSON' not in self._tool_inputs:
            raise InvalidToolInputError('Ref. selection JSON input is required (JSON)')
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('Ref. selection FASTA input is required (FASTA)')
        if 'DB' not in self._tool_inputs:
            raise InvalidToolInputError('Database input (DB) is required')
        super()._check_input()

    def _add_table_ref_metadata(self, section: HtmlReportSection, dir_db: Path, ref_by_seg: dict[str, pd.Series]) -> \
            None:
        """
        Adds a table with metadata for the selected reference.
        :param dir_db: Database directory
        :param ref_by_seg: Selected reference by segment
        :return: reference metadata dictionary
        """
        # Parse metadata file
        data_meta = dbutils.read_metadata(dir_db)

        # Create output table
        records_out = []
        for seg, ref in ref_by_seg.items():
            if ref is None:
                continue
            accession = seqidutils.get_accession(str(ref['ref_id']))
            metadata = data_meta[data_meta['id'] == accession].iloc[0]
            records_out.append({
                'Segment': seg,
                'Accession': accession,
                **{k.replace('_', ' ').lower().capitalize(): v
                   for k, v in metadata.fillna('n/a').items() if k not in ('id', 'accession')}
            })
        data_out = pd.DataFrame(records_out)

        # Add to report
        section.add_header('Reference metadata', 3)
        column_names = data_out.columns
        table_data = list(data_out.itertuples(index=False, name=None))
        # noinspection PyTypeChecker
        section.add_table(table_data, column_names, [('class', 'data')])

    def __add_overview_table(self, section: HtmlReportSection, ref_by_segment: dict[str, pd.Series]) -> None:
        """
        Adds a table with an overview.
        :param section: Report section
        :param ref_by_segment: Selected reference by segment
        :return: None
        """
        section.add_table([
            ['Database:', self._tool_inputs['DB'][0].path.name]
        ], table_attributes=[('class', 'information')])
        section.add_header('Overview', 3)
        table_data = [[
            seg,
            HtmlTableCell('Yes', color='green') if data is not None else HtmlTableCell('No', color='red'),
            data['ref_id_fmt'] if data is not None else '-',
            f"{int(str(data['median_mult'])):,}" if data is not None else '-',
            data['hashes'] if data is not None else '-'] for seg, data in ref_by_segment.items()
        ]
        header = ['Segment', 'Present', 'Ref. genome', 'Estimated cov.', 'Matching hashes']
        section.add_table(table_data, header, [('class', 'data')])

        # Add download link
        path_rel = Path('ref_selection', 'selected_ref_mash.fasta')
        section.add_file(self._tool_inputs['FASTA'][0].path, path_rel)
        section.add_link_to_file('Download reference genome (FASTA)', path_rel)

    def _add_segments_tables(self, section: HtmlReportSection, mash_combined: pd.DataFrame) -> None:
        """
        Adds a table with the results for each of the segments.
        :param section: HTML report section
        :param mash_combined: mash output combined across segments
        :return: None
        """
        section.add_header('Mash output (by segment)', 3)
        div_ = HtmlExpandableDiv('match_by_segment', 'Matches by segment')
        for seg, mash_out in mash_combined.groupby(by='segment'):
            div_.add_header(f'Segment {seg}', 3)
            column_keys = [c['key'] for c in ReporterRefSelection.COLS_SEGMENT]
            mash_out['identity'] = mash_out['identity'].apply(lambda x: f'{100 * x:.2f}')
            table_data = mash_out[column_keys].head(10)
            div_.add_table(
                HtmlTableFormatter.format_table_data(table_data, ReporterRefSelection.COLS_SEGMENT),
                [c['title'] for c in ReporterRefSelection.COLS_SEGMENT],
                [('class', 'data')])
        section.add_html_object(div_)
