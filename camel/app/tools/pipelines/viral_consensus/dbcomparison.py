from pathlib import Path

import pandas as pd
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import fastautils

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool
from camel.app.loggers import logger
from camel.app.toolkits.minimap2 import pafutils
from camel.app.toolkits.viralconsensus import dbutils


class DbComparison(Tool):
    """
    Compares the consensus sequence of each segment to the sequences in the reference database.
    The minimap2 PAF output is parsed and used to calculate the coverage and number of differences (SNPs, indels, etc).
    Metadat for the best matching sequences is also reported.
    """

    COLS_MATCHES = [
        'segment', 'rank', 'accession', 'coverage', 'coverage_called', 'identity', 'differences', 'mismatches',
        'indels', 'nb_accessions']
    COLS_EQUAL = ['segment', 'accession', 'unique_id', 'differences', 'identity']

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('viral consensus: database comparison', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        toolutils.check_input(self, keys_required=['DB'], keys_allowed=['PAF', 'FASTA', 'DB'])
        if 'PAF' in self._tool_inputs and 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('Consensus sequence input (FASTA) is required to process the alignments')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        dir_db = self._tool_inputs['DB'][0].path
        alignments = [DbComparison.parse_paf(io.path) for io in self._tool_inputs.get('PAF', [])]
        alignments = [data_aln for data_aln in alignments if len(data_aln) > 0]
        if len(alignments) > 0:
            data_bases = DbComparison._count_bases([io.path for io in self._tool_inputs['FASTA']])
            data_aln = DbComparison._combine_alignments(pd.concat(alignments, ignore_index=True), data_bases)
            matches, equal = [], []
            for segment, data_seg in data_aln.groupby('segment'):
                data_ranked = self._rank_matches(data_seg)
                data_clusters = DbComparison._parse_clusters(dir_db, segment)
                matches.append(DbComparison._add_nb_accessions(data_ranked, data_clusters))
                equal.append(DbComparison._expand_equally_close(data_ranked[data_ranked['rank'] == 1], data_clusters))
            data_matches = pd.concat(matches, ignore_index=True)
            data_equal = pd.concat(equal, ignore_index=True)
        else:
            logger.warning('No alignments to database sequences found')
            data_matches = pd.DataFrame(columns=DbComparison.COLS_MATCHES)
            data_equal = pd.DataFrame(columns=DbComparison.COLS_EQUAL)

        # Add the metadata
        data_meta = DbComparison._parse_metadata(dir_db, set(data_matches['accession']) | set(data_equal['accession']))
        data_matches = DbComparison._add_metadata(data_matches, data_meta)
        data_equal = DbComparison._add_metadata(data_equal, data_meta)

        # Create the output
        for key, filename, data in (
                ('TSV', 'closest_matches.tsv', data_matches),
                ('TSV_equal', 'equally_close_sequences.tsv', data_equal)):
            path_out = self.folder / filename
            data.to_csv(path_out, sep='\t', index=False)
            self._tool_outputs[key] = [ToolIOFile(path_out)]
        self._informs['min_coverage'] = float(self.get_param_value('min_coverage'))
        logger.info(
            f'Database comparison: {len(data_matches):,} matches reported, {len(data_equal):,} equally close '
            f'accessions')

    @staticmethod
    def parse_paf(path_paf: Path) -> pd.DataFrame:
        """
        Parses the minimap2 alignments of the database sequences (queries) against the consensus (target). Only the
        primary alignments are kept.
        :param path_paf: PAF file (with long cs tag)
        :return: Alignment statistics (one row per alignment)
        """
        data_paf = pafutils.read_paf(path_paf, tags=['cs'])
        data_out = pd.DataFrame({
            'segment': data_paf['target'].str.split('-').str[-1],
            'accession': data_paf['query'].str.rsplit('-', n=1).str[0]
        })
        data_stats = pd.DataFrame(
            [pafutils.count_differences(cs) for cs in data_paf['cs']], index=data_paf.index,
            columns=['mismatches', 'indels', 'compared_positions'])
        return data_out.join(data_stats)

    @staticmethod
    def _count_bases(paths_fasta: list[Path]) -> pd.DataFrame:
        """
        Counts the total and called (non-N) bases of the consensus sequence of each segment.
        :param paths_fasta: Consensus FASTA files
        :return: Total ('length') and called ('called') bases, indexed by segment
        """
        records = []
        for path_fasta in paths_fasta:
            for seq_id, seq in fastautils.read_as_dict(path_fasta).items():
                records.append({
                    'segment': seq_id.split('-')[-1],
                    'length': len(seq),
                    'called': len(seq) - str(seq.seq).upper().count('N')
                })
        return pd.DataFrame(records, columns=['segment', 'length', 'called']).set_index('segment')

    @staticmethod
    def _combine_alignments(data_aln: pd.DataFrame, data_bases: pd.DataFrame) -> pd.DataFrame:
        """
        Combines multiple alignments for each database sequence and calculates coverage metrics.

        Computes two coverage values per sequence:
        - `coverage`: Fraction of the total consensus sequence covered.
        - `coverage_called`: Fraction of called bases covered (excluding masked gap Ns).

        :param data_aln: Alignment statistics with one row per alignment.
        :param data_bases: Total and called base counts of the consensus, indexed by segment.
        :return: Aggregated alignment statistics with one row per database sequence.
        """
        data_out = data_aln.groupby(['segment', 'accession'], as_index=False, sort=False)[
            ['mismatches', 'indels', 'compared_positions']].sum()
        segments_missing = set(data_out['segment']) - set(data_bases.index)
        if len(segments_missing) > 0:
            raise ValueError(f"No consensus sequence for segment(s): {', '.join(sorted(segments_missing))}")
        compared = data_out['compared_positions']
        data_out['coverage'] = (compared / data_out['segment'].map(data_bases['length'])).clip(upper=1.0)
        data_out['coverage_called'] = (compared / data_out['segment'].map(data_bases['called'])).clip(upper=1.0)
        data_out['identity'] = (100 * (compared - data_out['mismatches']) / compared).where(compared > 0, 0.0)
        data_out['differences'] = data_out['mismatches'] + data_out['indels']
        return data_out[[
            'segment', 'accession', 'coverage', 'coverage_called', 'identity', 'differences', 'mismatches', 'indels']]

    def _rank_matches(self, data_seg: pd.DataFrame) -> pd.DataFrame:
        """
        Ranks the database sequences of a segment by the number of differences (ties share the rank). Sequences below
        the minimum coverage of the called bases are discarded, unless none of the sequences reaches it.
        :param data_seg: Alignment statistics for a single segment
        :return: Top matches (including all sequences tied with the closest match)
        """
        data_pass = data_seg[data_seg['coverage_called'] >= float(self.get_param_value('min_coverage'))]
        if len(data_pass) == 0:
            logger.warning(f"No matches above the minimum coverage ({data_seg['segment'].iloc[0]})")
            data_pass = data_seg
        data_ranked = data_pass.sort_values(
            by=['differences', 'coverage_called', 'identity'], ascending=[True, False, False], ignore_index=True)
        data_ranked.insert(1, 'rank', data_ranked['differences'].rank(method='min').astype(int))
        top_n = int(self.get_param_value('top_matches_reported'))
        return data_ranked[(data_ranked.index < top_n) | (data_ranked['rank'] == 1)]

    @staticmethod
    def _parse_clusters(dir_db: Path, segment: str) -> pd.DataFrame | None:
        """
        Parses the cluster table linking all accessions ('id') to the unique sequences in the database ('unique_id').
        :param dir_db: Database directory
        :param segment: Segment name
        :return: Cluster table, None if not available (databases created before v3)
        """
        path_clusters = dbutils.get_db_paths(dir_db, segment)['clusters']
        if path_clusters is None:
            return None
        return pd.read_table(path_clusters, usecols=['id', 'unique_id'], dtype=str)

    @staticmethod
    def _add_nb_accessions(data_ranked: pd.DataFrame, data_clusters: pd.DataFrame | None) -> pd.DataFrame:
        """
        Adds the number of accessions merged into each unique database sequence (missing without cluster table).
        :param data_ranked: Ranked matches of a single segment
        :param data_clusters: Cluster table (None if not available)
        :return: Ranked matches with number of accessions
        """
        if data_clusters is None:
            return data_ranked.assign(nb_accessions=pd.NA)
        nb_by_unique_id = data_clusters['unique_id'].value_counts()
        return data_ranked.assign(nb_accessions=data_ranked['accession'].map(nb_by_unique_id).fillna(1).astype(int))

    @staticmethod
    def _expand_equally_close(data_equal: pd.DataFrame, data_clusters: pd.DataFrame | None) -> pd.DataFrame:
        """
        Expands the equally close unique database sequences to all accessions merged into them. Without cluster table,
        only the unique sequences are reported.
        :param data_equal: Equally close matches (rank 1) of a single segment
        :param data_clusters: Cluster table (None if not available)
        :return: Equally close accessions
        """
        data_out = data_equal[['segment', 'accession', 'differences', 'identity']].rename(
            columns={'accession': 'unique_id'})
        if data_clusters is None:
            data_out['accession'] = data_out['unique_id']
        else:
            data_out = data_out.merge(data_clusters.rename(columns={'id': 'accession'}), on='unique_id', how='left')
            data_out['accession'] = data_out['accession'].fillna(data_out['unique_id'])
        return data_out[DbComparison.COLS_EQUAL]

    @staticmethod
    def _parse_metadata(dir_db: Path, accessions: set[str]) -> pd.DataFrame:
        """
        Parses the metadata of the given accessions.
        :param dir_db: Database directory
        :param accessions: Accessions
        :return: Metadata
        """
        data_meta = pd.read_table(
            dir_db / 'sequence_metadata.tsv', dtype=str, keep_default_na=False, na_values=['-', ''])
        return data_meta[data_meta['id'].isin(accessions)]

    @staticmethod
    def _add_metadata(data: pd.DataFrame, data_meta: pd.DataFrame) -> pd.DataFrame:
        """
        Adds the metadata of the database sequences, based on the 'accession' column.
        :param data: Input data
        :param data_meta: Metadata
        :return: Data with metadata
        """
        return data.merge(data_meta, left_on='accession', right_on='id', how='left').drop(columns='id')
