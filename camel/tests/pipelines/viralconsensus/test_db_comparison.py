import unittest
from pathlib import Path

import pandas as pd
from camelcore.app.io.tooliodirectory import ToolIODirectory
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.pipelines.viral_consensus.dbcomparison import DbComparison
from camel.app.tools.pipelines.viral_consensus.reporterdbcomparison import ReporterDbComparison


class TestDbComparison(CamelTestSuite):
    """
    Tests for the database comparison. The test database contains segments HA and PB2 (no result), ACC1b is identical
    to ACC1 (v3 database with cluster table and manifest, v2 database without).
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines', 'viral_consensus', 'db_comparison')
    DIR_DB = test_file_dir / 'db_v3'
    DIR_DB_V2 = test_file_dir / 'db_v2'
    FASTA = test_file_dir / 'consensus.fasta'
    FASTA_GAP = test_file_dir / 'consensus_gap.fasta'
    PAF_HA = test_file_dir / 'HA.paf'
    PAF_HA_GAP = test_file_dir / 'HA-gap.paf'
    PAF_NA = test_file_dir / 'NA.paf'
    PAF_EMPTY = test_file_dir / 'empty.paf'

    def _run_tool(
            self, path_paf: Path = PAF_HA, dir_db: Path = DIR_DB, path_fasta: Path = FASTA, **params) -> DbComparison:
        """
        Runs the database comparison.
        :param path_paf: PAF file
        :param dir_db: Database directory
        :param path_fasta: Consensus FASTA file
        :param params: Tool parameters
        :return: Tool
        """
        db_comparison = DbComparison()
        db_comparison.add_input_files({
            'PAF': [ToolIOFile(path_paf)],
            'FASTA': [ToolIOFile(path_fasta)],
            'DB': [ToolIODirectory(dir_db)]
        })
        if params:
            db_comparison.update_parameters(**params)
        db_comparison.run(self.running_dir)
        self.verify_output_files(db_comparison, 'TSV')
        self.verify_output_files(db_comparison, 'TSV_equal')
        return db_comparison

    def _run_reporter(self, db_comparison: DbComparison, dir_db: Path = DIR_DB, **params) -> str:
        """
        Creates the report section for the database comparison output.
        :param db_comparison: Database comparison tool
        :param dir_db: Database directory
        :param params: Reporter parameters
        :return: HTML of the report section
        """
        reporter = ReporterDbComparison()
        reporter.add_input_files({
            'TSV': db_comparison.tool_outputs['TSV'],
            'TSV_equal': db_comparison.tool_outputs['TSV_equal'],
            'DB': [ToolIODirectory(dir_db)]
        })
        reporter.add_input_informs({'db_comparison': db_comparison.informs})
        if params:
            reporter.update_parameters(**params)
        reporter.run(self.running_dir)
        output_section = reporter.tool_outputs['VAL_HTML'][0].value
        CamelTestSuite.export_report_section(output_section, self.running_dir / 'report')
        return output_section.to_html()

    def test_db_comparison(self) -> None:
        """
        Tests the database comparison.
        :return: None
        """
        db_comparison = self._run_tool()

        # ACC1 and ACC2 are identical (ACC1 has the higher coverage), ACC4 is below the minimum coverage
        data_matches = pd.read_table(db_comparison.tool_outputs['TSV'][0].path)
        self.assertEqual(list(data_matches['accession']), ['ACC1', 'ACC2', 'ACC3'])
        self.assertEqual(list(data_matches['rank']), [1, 1, 3])
        self.assertEqual(list(data_matches['nb_accessions']), [2, 1, 1])
        self.assertEqual(list(data_matches['Host']), ['Homo sapiens', 'Sus scrofa', 'Homo sapiens'])

        # Equally close sequences expanded to all accessions (ACC1b is identical to ACC1)
        data_equal = pd.read_table(db_comparison.tool_outputs['TSV_equal'][0].path)
        self.assertEqual(list(data_equal['accession']), ['ACC1', 'ACC1b', 'ACC2'])
        self.assertEqual(list(data_equal['unique_id']), ['ACC1', 'ACC1', 'ACC2'])
        self.assertEqual(list(data_equal['Host']), ['Homo sapiens', 'Gallus gallus', 'Sus scrofa'])

    def test_db_comparison_split_alignment(self) -> None:
        """
        Tests a consensus with a masked gap (N stretch): the alignments on both sides of the gap are combined and the
        coverage is calculated on the called bases.
        :return: None
        """
        db_comparison = self._run_tool(path_paf=TestDbComparison.PAF_HA_GAP, path_fasta=TestDbComparison.FASTA_GAP)
        data_matches = pd.read_table(db_comparison.tool_outputs['TSV'][0].path)
        self.assertEqual(list(data_matches['accession']), ['ACC1'])
        self.assertAlmostEqual(data_matches['coverage'].iloc[0], 40 / 60)
        self.assertEqual(data_matches['coverage_called'].iloc[0], 1.0)
        self.assertEqual(data_matches['differences'].iloc[0], 1)
        self.assertEqual(data_matches['identity'].iloc[0], 97.5)

        # The consensus with N positions is flagged in the report
        html = self._run_reporter(db_comparison)
        self.assertIn('<td class="red">66.67%</td>', html)
        self.assertIn('The closest match covers less than 90% of the consensus sequence for segment(s) HA', html)

    def test_db_comparison_no_clusters(self) -> None:
        """
        Tests the database comparison with a database without cluster table (created before v3).
        :return: None
        """
        db_comparison = self._run_tool(dir_db=TestDbComparison.DIR_DB_V2)
        data_matches = pd.read_table(db_comparison.tool_outputs['TSV'][0].path)
        self.assertTrue(data_matches['nb_accessions'].isna().all())
        data_equal = pd.read_table(db_comparison.tool_outputs['TSV_equal'][0].path)
        self.assertEqual(list(data_equal['accession']), ['ACC1', 'ACC2'])

    def test_db_comparison_low_coverage(self) -> None:
        """
        Tests that the best match is reported when no match reaches the minimum coverage.
        :return: None
        """
        db_comparison = self._run_tool(min_coverage=1.01)
        data_matches = pd.read_table(db_comparison.tool_outputs['TSV'][0].path)
        self.assertEqual(data_matches['accession'].iloc[0], 'ACC1')

    def test_db_comparison_segment_na(self) -> None:
        """
        Tests that the NA segment is not parsed as a missing value.
        :return: None
        """
        html = self._run_reporter(self._run_tool(path_paf=TestDbComparison.PAF_NA))
        self.assertIn('<td>NA</td><td>ACC1</td>', html)
        self.assertIn('Segment NA', html)

    def test_db_comparison_no_alignments(self) -> None:
        """
        Tests the database comparison and report when there are no alignments (empty PAF file).
        :return: None
        """
        db_comparison = self._run_tool(path_paf=TestDbComparison.PAF_EMPTY)
        self.assertEqual(len(pd.read_table(db_comparison.tool_outputs['TSV'][0].path)), 0)
        self.assertEqual(len(pd.read_table(db_comparison.tool_outputs['TSV_equal'][0].path)), 0)
        self.assertIn('No matching database sequences found', self._run_reporter(db_comparison))

    def test_reporter(self) -> None:
        """
        Tests the database comparison reporter.
        :return: None
        """
        html = self._run_reporter(self._run_tool())
        # Database information from the manifest
        self.assertIn('<td>Database:</td><td>influenza_a-test</td>', html)
        self.assertIn('<td>Download date:</td><td>2026-09-28</td>', html)
        # Overview: 2 equally close sequences, 3 accessions, metadata summarized for the 3 accessions
        self.assertIn(
            '<td class="green">100.00%</td><td class="green">100.00%</td><td>0</td><td>2</td><td>3</td>'
            '<td>H3N2 (2), H1N1 (1)</td>'
            '<td>Gallus gallus (1), Homo sapiens (1), Sus scrofa (1)</td>'
            '<td>Belgium (1), France (1), n/a (1)</td><td>2019 - 2021-03</td>', html)
        # Metadata and number of accessions of all top matches, missing values shown as 'n/a'
        self.assertIn('<td>ACC2</td><td>100.00%</td><td>97.50%</td><td>0</td><td>0</td><td>0</td><td>1</td>'
                      '<td>H1N1</td><td>Sus scrofa</td><td>n/a</td><td>2019</td>', html)
        self.assertIn('equally_close_sequences.tsv', html)
        # PB2 is in the database but has no result: row with 0.00% aligned and alert
        self.assertIn(
            '<td>PB2</td><td>n/a</td><td class="red">n/a</td><td class="red">0.00%</td>'
            '<td>n/a</td><td>n/a</td><td>n/a</td><td>n/a</td><td>n/a</td><td>n/a</td><td>n/a</td>', html)
        self.assertIn('No database comparison for segment(s) PB2', html)
        # The closest HA match covers the full consensus and is closely related
        self.assertNotIn('The closest match covers less than', html)
        self.assertNotIn('No closely related database sequence was found', html)

    def test_reporter_low_coverage(self) -> None:
        """
        Tests that an alert is shown when the closest match does not cover the consensus sequence.
        :return: None
        """
        html = self._run_reporter(self._run_tool(), coverage_warn=102.0, coverage_fail=101.0)
        self.assertIn('<td class="red">100.00%</td>', html)
        self.assertIn('The closest match covers less than 101% of the consensus sequence for segment(s) HA', html)

    def test_reporter_summarize_metadata(self) -> None:
        """
        Tests the summaries of the metadata of the equally close accessions.
        :return: None
        """
        values = pd.Series(['USA'] * 10 + ['New Zealand'] * 4 + ['Canada'] * 2 + ['Belgium', 'France', None])
        self.assertEqual(
            ReporterDbComparison._summarize_counts(values), 'USA (10), New Zealand (4), Canada (2), other (3)')
        self.assertEqual(ReporterDbComparison._summarize_dates(pd.Series(['2019-02-05', '2018-10', None])),
                         '2018-10 - 2019-02')
        self.assertEqual(ReporterDbComparison._summarize_dates(pd.Series([None], dtype=str)), 'n/a')

    def test_reporter_no_clusters(self) -> None:
        """
        Tests the report for a database without cluster table and manifest (created before v3).
        :return: None
        """
        dir_db = TestDbComparison.DIR_DB_V2
        html = self._run_reporter(self._run_tool(dir_db=dir_db), dir_db=dir_db)
        self.assertIn('<td class="green">100.00%</td><td>0</td><td>2</td><td>n/a</td>', html)
        # No manifest: database directory name, download date and number of sequences not available
        self.assertIn('<td>Database:</td><td>db_v2</td>', html)
        self.assertIn('<td>Download date:</td><td>n/a</td>', html)
        self.assertIn('<td>Nb. sequences:</td><td>n/a</td>', html)

    def test_reporter_low_identity(self) -> None:
        """
        Tests that an alert is shown when the closest match is below the identity threshold.
        :return: None
        """
        html = self._run_reporter(self._run_tool(), identity_warn=101.0, identity_fail=100.5)
        self.assertIn('No closely related database sequence was found', html)


if __name__ == '__main__':
    unittest.main()
