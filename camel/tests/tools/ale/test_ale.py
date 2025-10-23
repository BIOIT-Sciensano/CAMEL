import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.ale.ale import ALE
from camel.app.tools.ale.ale2wiggle import ALE2Wiggle


class TestALE(CamelTestSuite):
    """
    Initializes the ALE testing tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('ale')
    FILE_SAM = ToolIOFile(test_file_dir / 'S1_ilmn_ppolish.sam')
    FILE_FASTA = ToolIOFile(test_file_dir / 'sequences.fasta')
    FILE_ALE = ToolIOFile(test_file_dir / 'ALE_ilmn.ale')

    def test_ale(self) -> None:
        """
        Actually testing ALE from illumina sequencing data
        """
        ale = ALE()
        ale.add_input_files({'FASTA': [TestALE.FILE_FASTA], 'SAM': [TestALE.FILE_SAM]})
        ale.run(self.running_dir)
        self.verify_output_files(ale, 'ALE')

    def test_ale2wiggle(self) -> None:
        """
        Actually testing ALE2Wiggle on ALE output file
        """
        ale2wiggle = ALE2Wiggle()
        ale2wiggle.add_input_files({'ALE': [TestALE.FILE_ALE]})
        ale2wiggle.run(self.running_dir)
        self.verify_output_files(ale2wiggle, 'TSV', 4)


if __name__ == '__main__':
    unittest.main()
