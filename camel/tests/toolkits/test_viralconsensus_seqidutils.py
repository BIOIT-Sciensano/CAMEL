import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.toolkits.viralconsensus import seqidutils


class TestViralConsensusSeqIdUtils(CamelTestSuite):
    """
    Tests the viral consensus sequence id utilities.
    """

    def test_get_segment(self) -> None:
        """
        Tests extracting the segment from a sequence id.
        :return: None
        """
        self.assertEqual(seqidutils.get_segment('OQ839344.1-HA'), 'HA')
        self.assertEqual(seqidutils.get_segment('sample-2024-01-PB2'), 'PB2')

    def test_get_accession(self) -> None:
        """
        Tests extracting the accession from a sequence id.
        :return: None
        """
        self.assertEqual(seqidutils.get_accession('OQ839344.1-HA'), 'OQ839344.1')
        self.assertEqual(seqidutils.get_accession('sample-2024-01-PB2'), 'sample-2024-01')
        self.assertEqual(seqidutils.get_accession('OQ839344.1'), 'OQ839344.1')

    def test_is_valid_seq_id(self) -> None:
        """
        Tests the validation of the sequence id format.
        :return: None
        """
        self.assertTrue(seqidutils.is_valid_seq_id('OQ839344.1-HA'))
        self.assertFalse(seqidutils.is_valid_seq_id('OQ839344.1'))


if __name__ == '__main__':
    unittest.main()
