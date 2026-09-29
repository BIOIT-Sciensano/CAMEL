import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.toolkits.minimap2 import pafutils


class TestPafUtils(CamelTestSuite):
    """
    Tests the pafutils module.
    """

    PAF_RECORDS = [
        ['seq1', 40, 0, 40, '+', 'target', 40, 0, 40, 40, 40, 60, 'NM:i:0', 'tp:A:P', 'cs:Z:=ACGT'],
        ['seq2', 40, 0, 40, '-', 'target', 40, 0, 40, 39, 40, 60, 'tp:A:S', 'cs:Z:=AC*ag=T'],
        ['seq3', 40, 0, 40, '+', 'target', 40, 0, 40, 39, 40, 60, 'tp:A:P'],
    ]

    def _write_paf(self) -> None:
        """
        Writes the test PAF file.
        :return: None
        """
        self.path_paf = self.running_dir / 'test.paf'
        self.path_paf.write_text(''.join('\t'.join(str(x) for x in r) + '\n' for r in TestPafUtils.PAF_RECORDS))

    def test_read_paf(self) -> None:
        """
        Tests reading a PAF file (primary alignments only).
        :return: None
        """
        self._write_paf()
        data_paf = pafutils.read_paf(self.path_paf, tags=['cs', 'NM'])
        self.assertEqual(list(data_paf['query']), ['seq1', 'seq3'])
        self.assertEqual(data_paf['cs'].iloc[0], '=ACGT')
        self.assertEqual(data_paf['NM'].iloc[0], '0')
        self.assertTrue(data_paf[['cs', 'NM']].iloc[1].isna().all())
        self.assertEqual(list(data_paf['target_end']), [40, 40])

    def test_read_paf_all_alignments(self) -> None:
        """
        Tests reading a PAF file, including the secondary alignments.
        :return: None
        """
        self._write_paf()
        data_paf = pafutils.read_paf(self.path_paf, tags=['cs'], primary_only=False)
        self.assertEqual(list(data_paf['query']), ['seq1', 'seq2', 'seq3'])

    def test_read_paf_empty(self) -> None:
        """
        Tests reading an empty PAF file.
        :return: None
        """
        path_paf = self.running_dir / 'empty.paf'
        path_paf.touch()
        data_paf = pafutils.read_paf(path_paf, tags=['cs'])
        self.assertEqual(len(data_paf), 0)
        self.assertIn('cs', data_paf.columns)

    def test_count_differences(self) -> None:
        """
        Tests counting the differences from the cs tag.
        :return: None
        """
        # Identical
        stats = pafutils.count_differences('=ACGTACGT')
        self.assertEqual((stats['differences'], stats['identity']), (0, 100))
        # Mismatch
        stats = pafutils.count_differences('=ACG*ta=ACGT')
        self.assertEqual((stats['mismatches'], stats['differences']), (1, 1))
        # N / ambiguous bases (both reported as 'n' by minimap2) are not compared
        stats = pafutils.count_differences('=ACG*na*an=GT')
        self.assertEqual((stats['differences'], stats['identity']), (0, 100))
        # A 3 bp deletion / insertion is a single difference
        self.assertEqual(pafutils.count_differences('=ACG-ttt=ACG')['differences'], 1)
        self.assertEqual(pafutils.count_differences('=ACG+ttt=ACG')['differences'], 1)


if __name__ == '__main__':
    unittest.main()
