import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.unicycler.unicycler import Unicycler
from camel.tests import longRunningTest


class TestUnicycler(CamelTestSuite):
    """
    Contains tests for the Unicycler tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('unicycler')
    FILE_FQ = ToolIOFile(test_file_dir / 'subsampled_ont.fastq.gz')
    FILE_FQ_1 = ToolIOFile(test_file_dir / 'pilsner_subsampled_1.fastq.gz')
    FILE_FQ_2 = ToolIOFile(test_file_dir / 'pilsner_subsampled_2.fastq.gz')

    def test_unicycler_check_input(self) -> None:
        """
        Tests the Unicycler check_input method.
        :return: None
        """
        unicycler = Unicycler()

        # Emtpy input
        unicycler.add_input_files({})
        with self.assertRaises(InvalidToolInputError):
            unicycler._check_input()

        # Invalid input type
        unicycler.add_input_files({'BAM': [TestUnicycler.FILE_FQ]})
        with self.assertRaises(InvalidToolInputError):
            unicycler._check_input()

        # Valid input
        unicycler.add_input_files({'FASTQ_SE': [TestUnicycler.FILE_FQ]})
        unicycler._check_input()

    @longRunningTest()
    def test_unicycler(self) -> None:
        """
        Tests Unicycler 0.5.0.
        """
        unicycler = Unicycler()
        unicycler.add_input_files({
            'FASTQ_SE': [TestUnicycler.FILE_FQ],
            'FASTQ_PE': [TestUnicycler.FILE_FQ_1, TestUnicycler.FILE_FQ_2]})
        unicycler.run(self.running_dir)
        self.verify_output_files(unicycler, 'FASTA')


if __name__ == '__main__':
    unittest.main()
