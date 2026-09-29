from pathlib import Path

import pandas as pd
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.tool import Tool
from camel.app.loggers import logger
from camel.app.toolkits.minimap2 import pafutils


class DbComparison(Tool):
    """
    Compares the consensus sequence of each segment to the sequences in the reference database. The alignments of the
    database sequences (minimap2 PAF output, consensus as target, with long cs tag) are ranked by the number of
    differences and linked to the database metadata.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('viral consensus: database comparison', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        toolutils.check_input(self, keys_required=['PAF', 'DB'])
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        data_aln = pd.concat([DbComparison.parse_paf(io.path) for io in self._tool_inputs['PAF']], ignore_index=True)
        data_matches = pd.concat(
            [self._rank_matches(data_seg) for _, data_seg in data_aln.groupby('segment')], ignore_index=True)
        data_matches = DbComparison._add_metadata(data_matches, self._tool_inputs['DB'][0].path)

        path_out = self.folder / 'closest_matches.tsv'
        data_matches.to_csv(path_out, sep='\t', index=False, na_rep='NA')
        self._tool_outputs['TSV'] = [ToolIOFile(path_out)]
        self._informs['min_coverage'] = float(self.get_param_value('min_coverage'))
        logger.info(f'Database comparison: {len(data_matches):,} matches reported')

    @staticmethod
    def parse_paf(path_paf: Path) -> pd.DataFrame:
        """
        Parses the minimap2 alignments of the database sequences (queries) against the consensus (target). Only the
        primary alignments are kept.
        :param path_paf: PAF file (with long cs tag)
        :return: Alignment statistics
        """
        data_paf = pafutils.read_paf(path_paf, tags=['cs'])
        data_out = pd.DataFrame({
            'segment': data_paf['target'].str.split('-').str[-1],
            'accession': data_paf['query'].str.rsplit('-', n=1).str[0],
            'coverage': (data_paf['target_end'] - data_paf['target_start']) / data_paf['target_length']
        })
        data_stats = pd.DataFrame([pafutils.count_differences(cs) for cs in data_paf['cs']], index=data_paf.index)
        return data_out.join(data_stats)

    def _rank_matches(self, data_seg: pd.DataFrame) -> pd.DataFrame:
        """
        Ranks the database sequences of a segment by the number of differences (ties share the rank). Sequences below
        the minimum coverage are discarded, unless none of the sequences reaches it.
        :param data_seg: Alignment statistics for a single segment
        :return: Top matches (including all sequences tied with the closest match)
        """
        data_pass = data_seg[data_seg['coverage'] >= float(self.get_param_value('min_coverage'))]
        if len(data_pass) == 0:
            logger.warning(f"No matches above the minimum coverage ({data_seg['segment'].iloc[0]})")
            data_pass = data_seg
        data_ranked = data_pass.sort_values(
            by=['differences', 'coverage', 'identity'], ascending=[True, False, False], ignore_index=True)
        data_ranked.insert(1, 'rank', data_ranked['differences'].rank(method='min').astype(int))
        top_n = int(self.get_param_value('top_matches_reported'))
        return data_ranked[(data_ranked.index < top_n) | (data_ranked['rank'] == 1)]

    @staticmethod
    def _add_metadata(data_matches: pd.DataFrame, dir_db: Path) -> pd.DataFrame:
        """
        Adds the metadata of the database sequences.
        :param data_matches: Matches
        :param dir_db: Database directory
        :return: Matches with metadata
        """
        data_meta = pd.read_table(
            dir_db / 'sequence_metadata.tsv', dtype=str, keep_default_na=False, na_values=['-', ''])
        data_meta = data_meta[data_meta['id'].isin(data_matches['accession'])]
        return data_matches.merge(data_meta, left_on='accession', right_on='id', how='left').drop(columns='id')
