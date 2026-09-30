import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.toolkits.viralconsensus import dbutils


class TestViralConsensusDbUtils(CamelTestSuite):
    """
    Tests the viral consensus database utilities (v3 database with manifest and cluster tables, v2 database without).
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines', 'viral_consensus', 'db_comparison')
    DIR_DB = test_file_dir / 'db_v3'
    DIR_DB_V2 = test_file_dir / 'db_v2'

    def test_get_db_paths(self) -> None:
        """
        Tests resolving the database paths (cluster table only available for v3 databases).
        :return: None
        """
        db_paths = dbutils.get_db_paths(TestViralConsensusDbUtils.DIR_DB, 'HA')
        self.assertEqual(db_paths['fasta'], TestViralConsensusDbUtils.DIR_DB / 'fasta_by_segment' / 'HA.fasta')
        self.assertEqual(db_paths['mash'], TestViralConsensusDbUtils.DIR_DB / 'mash' / 'HA.msh')
        self.assertEqual(db_paths['clusters'], TestViralConsensusDbUtils.DIR_DB / 'clusters' / 'HA.tsv')
        self.assertFalse(db_paths['is_full'])
        self.assertIsNone(dbutils.get_db_paths(TestViralConsensusDbUtils.DIR_DB_V2, 'HA')['clusters'])

    def test_get_segments(self) -> None:
        """
        Tests retrieving the segments in the database.
        :return: None
        """
        self.assertEqual(dbutils.get_segments(TestViralConsensusDbUtils.DIR_DB), ['HA', 'PB2'])

    def test_read_genome_info(self) -> None:
        """
        Tests reading the genome information.
        :return: None
        """
        self.assertEqual(dbutils.read_genome_info(TestViralConsensusDbUtils.DIR_DB)['genome_size'], 13500)

    def test_read_manifest(self) -> None:
        """
        Tests reading the manifest (empty for databases without manifest).
        :return: None
        """
        manifest = dbutils.read_manifest(TestViralConsensusDbUtils.DIR_DB)
        self.assertEqual(manifest['ncbi_download_date'], '2026-09-28')
        self.assertEqual(dbutils.read_manifest(TestViralConsensusDbUtils.DIR_DB_V2), {})

    def test_read_metadata(self) -> None:
        """
        Tests reading the sequence metadata ('-' parsed as missing value, dates read as text).
        :return: None
        """
        data_meta = dbutils.read_metadata(TestViralConsensusDbUtils.DIR_DB)
        self.assertEqual(len(data_meta), 4)
        data_meta = dbutils.read_metadata(TestViralConsensusDbUtils.DIR_DB, accessions={'ACC2', 'ACC3'})
        self.assertEqual(list(data_meta['id']), ['ACC2', 'ACC3'])
        self.assertTrue(data_meta['Country'].isna().iloc[0])
        self.assertEqual(data_meta['Collection_Date'].iloc[0], '2019')

    def test_read_seq_lengths(self) -> None:
        """
        Tests reading the sequence lengths from the FASTA index files.
        :return: None
        """
        self.assertEqual(dbutils.read_seq_lengths(TestViralConsensusDbUtils.DIR_DB), {'HA': {'ACC1-HA': 40}})


if __name__ == '__main__':
    unittest.main()
