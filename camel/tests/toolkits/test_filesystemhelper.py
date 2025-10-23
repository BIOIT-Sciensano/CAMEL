import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.utils import fileutils


class TestFileSystemHelper(CamelTestSuite):
    """
    Tests the Fastq utils module.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('components')

    def test_is_gzipped_uncompressed(self) -> None:
        """
        Tests the get sample name function for MiSEQ format.
        :return: None
        """
        path = TestFileSystemHelper.test_file_dir / 'fq-file.fastq'
        self.assertFalse(fileutils.is_gzipped(path))

    def test_is_gzipped_compressed(self) -> None:
        """
        Tests the get sample name function for MiSEQ format.
        :return: None
        """
        path = TestFileSystemHelper.test_file_dir / 'fq-file.fastq.gz'
        self.assertTrue(fileutils.is_gzipped(path))

    def test_is_gzipped_compressed_no_ext(self) -> None:
        """
        Tests the get sample name function for MiSEQ format.
        :return: None
        """
        path = TestFileSystemHelper.test_file_dir / 'fq-file-no-ext'
        self.assertFalse(fileutils.is_gzipped(path))

    def test_is_gzipped_uncompressed_no_ext(self) -> None:
        """
        Tests the get sample name function for MiSEQ format.
        :return: None
        """
        path = TestFileSystemHelper.test_file_dir / 'fq-file-no-ext-gzipped'
        self.assertTrue(fileutils.is_gzipped(path))


if __name__ == '__main__':
    unittest.main()
