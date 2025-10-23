import unittest

from camel.app.core.utils import sambamutils
from camel.app.core.cameltestsuite import CamelTestSuite


class TestBAMUtils(CamelTestSuite):
    """
    Tests the BAM utils module.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('components')

    def test_get_record_count_bam(self) -> None:
        """
        Tests the get_record_count function on BAM input.
        :return: None
        """
        nb_reads = sambamutils.get_record_count(TestBAMUtils.test_file_dir / 'toy.bam')
        self.assertGreater(nb_reads, 0)

    def test_get_record_count_sam(self) -> None:
        """
        Tests the get_record_count function on SAM input.
        :return: None
        """
        nb_reads = sambamutils.get_record_count(TestBAMUtils.test_file_dir / 'toy.sam')
        self.assertGreater(nb_reads, 0)

    def test_is_empty_count_sam(self) -> None:
        """
        Tests the is_empty function on SAM input.
        :return: None
        """
        self.assertFalse(sambamutils.is_empty(TestBAMUtils.test_file_dir / 'toy.sam'))

    def test_is_empty_count_bam(self) -> None:
        """
        Tests the is_empty function on SAM input.
        :return: None
        """
        self.assertFalse(sambamutils.is_empty(TestBAMUtils.test_file_dir / 'toy.bam'))

    def test_create_empty_sam(self) -> None:
        """
        Tests the create_empty function with SAM output.
        :return: None
        """
        path_out = self.running_dir / 'output.sam'
        path_ref = TestBAMUtils.test_file_dir / 'toy.fasta'
        sambamutils.create_empty(path_out=path_out, fasta_in=path_ref, compress=False)
        self.assertTrue(sambamutils.is_empty(path_out))

    def test_create_empty_bam(self) -> None:
        """
        Tests the create_empty function with SAM output.
        :return: None
        """
        path_out = self.running_dir / 'output.bam'
        path_ref = TestBAMUtils.test_file_dir / 'toy.fasta'
        sambamutils.create_empty(path_out=path_out, fasta_in=path_ref, compress=True)
        self.assertTrue(sambamutils.is_empty(path_out))


if __name__ == '__main__':
    unittest.main()
