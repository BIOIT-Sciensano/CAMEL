import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.prodigal.prodigal import Prodigal


class TestProdigal(CamelTestSuite):
    """
    Tests the Prodigal tool.
    """
    # Get test file and reference file directories
    test_file_dir = CamelTestSuite.get_test_file_dir('prodigal')
    input_fasta = test_file_dir / 'ecoli_contigs.fasta'

    def test_prodigal(self) -> None:
        """
        Tests the Prodigal tool.
        :return: None
        """
        prodigal = Prodigal()
        prodigal.add_input_files({'FASTA': [ToolIOFile(TestProdigal.input_fasta)]})
        prodigal.run(self.running_dir)
        self.verify_output_files(prodigal, 'GBK')
        self.verify_output_files(prodigal, 'FASTA')


if __name__ == '__main__':
    unittest.main()
