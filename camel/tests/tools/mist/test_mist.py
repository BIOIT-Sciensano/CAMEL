import unittest

from camelcore.app.io.tooliodirectory import ToolIODirectory
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mist.mistcall import MiSTCall
from camel.app.tools.mist.mistlincode import MiSTLincode


class TestMiST(CamelTestSuite):
    """
    Tests the MiST tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('typing')
    input_fasta = test_file_dir / 'neisseria_mc58.fasta'
    input_db = test_file_dir / 'scheme_mlst_neisseria'

    def _verify_output(self, tool_instance: MiSTCall) -> None:
        """
        Verifies that the tool output is correct.
        :param tool_instance: The tool instance.
        :return: None
        """
        for k in ('loci_detected', 'loci_total'):
            self.assertIn(k, tool_instance.informs)
        self.verify_output_files(tool_instance, 'JSON')

    def test_mist_call(self) -> None:
        """
        Tests the MiST call function with a MLST scheme
        :return: None
        """
        mist_call = MiSTCall()
        mist_call.add_input_files({
            'FASTA': [ToolIOFile(TestMiST.test_file_dir / 'neisseria_mc58.fasta')],
            'DB': [ToolIODirectory(TestMiST.input_db / 'mist')]
        })
        mist_call.run(self.running_dir)
        self._verify_output(mist_call)

    def test_mist_call_new_allele(self) -> None:
        """
        Tests the MiST call function with a MLST scheme
        :return: None
        """
        mist_call = MiSTCall()
        mist_call.add_input_files({
            'FASTA': [ToolIOFile(TestMiST.test_file_dir / 'neisseria_mc58_new_allele.fasta')],
            'DB': [ToolIODirectory(TestMiST.input_db / 'mist')]
        })
        mist_call.run(self.running_dir)
        self._verify_output(mist_call)

    def test_mist_lincode(self) -> None:
        """
        Tests MiST lincode.
        """
        mist_lincode = MiSTLincode()
        mist_lincode.add_input_files({
            'JSON': [ToolIOFile(CamelTestSuite.get_test_file_dir('mist') / 'mist_out_klebsiella_scgmlst.json')],
            'DB': [ToolIODirectory(CamelTestSuite.get_test_file_dir('mist') / 'dummy_db')],
        })
        mist_lincode.run(self.running_dir)
        self.verify_output_files(mist_lincode, 'JSON')
        for key in ('lincode_full', 'lincode_partial', 'lincode_partial_str', 'thresholds', ):
            self.assertIn(key, mist_lincode.informs)


if __name__ == '__main__':
    unittest.main()
