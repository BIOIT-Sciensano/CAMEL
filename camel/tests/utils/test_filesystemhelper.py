import unittest

from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.testing.cameltestsuite import CamelTestSuite


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
        self.assertFalse(FileSystemHelper.is_gzipped(path))

    def test_is_gzipped_compressed(self) -> None:
        """
        Tests the get sample name function for MiSEQ format.
        :return: None
        """
        path = TestFileSystemHelper.test_file_dir / 'fq-file.fastq.gz'
        self.assertTrue(FileSystemHelper.is_gzipped(path))

    def test_is_gzipped_compressed_no_ext(self) -> None:
        """
        Tests the get sample name function for MiSEQ format.
        :return: None
        """
        path = TestFileSystemHelper.test_file_dir / 'fq-file-no-ext'
        self.assertFalse(FileSystemHelper.is_gzipped(path))

    def test_is_gzipped_uncompressed_no_ext(self) -> None:
        """
        Tests the get sample name function for MiSEQ format.
        :return: None
        """
        path = TestFileSystemHelper.test_file_dir / 'fq-file-no-ext-gzipped'
        self.assertTrue(FileSystemHelper.is_gzipped(path))


if __name__ == '__main__':
    unittest.main()
