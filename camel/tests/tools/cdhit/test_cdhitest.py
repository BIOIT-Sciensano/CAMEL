import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.cdhit.cdhitest import CDHitEst


class TestCDHitEst(CamelTestSuite):
    """
    Tests the CD-HIT-EST tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('cdhit')
    FILE_FASTA = ToolIOFile(test_file_dir / 'seqs_in.fasta')

    def test_cdhitest(self) -> None:
        """
        Tests CD-HIT-EST with default parameters.
        """
        tool = CDHitEst()
        tool.add_input_files({'FASTA': [TestCDHitEst.FILE_FASTA]})
        tool.run(self.running_dir)
        self.verify_output_files(tool, 'FASTA')
        self.assertIn('clusters', tool.informs)
        self.assertGreater(len(tool.informs['clusters']), 0)

    def test_cdhitest_identity_threshold(self) -> None:
        """
        Tests CD-HIT-EST with a lower identity threshold, expecting more sequences to be clustered.
        """
        tool_high = CDHitEst()
        tool_high.add_input_files({'FASTA': [TestCDHitEst.FILE_FASTA]})
        tool_high.update_parameters(identity_threshold='0.9999')
        tool_high.run(self.running_dir)

        tool_low = CDHitEst()
        tool_low.add_input_files({'FASTA': [TestCDHitEst.FILE_FASTA]})
        tool_low.update_parameters(identity_threshold='0.8')
        tool_low.run(self.running_dir)

        self.assertGreaterEqual(
            len(tool_high.informs['clusters']),
            len(tool_low.informs['clusters']),
            "Lower identity threshold should produce fewer or equal clusters"
        )


if __name__ == '__main__':
    unittest.main()
