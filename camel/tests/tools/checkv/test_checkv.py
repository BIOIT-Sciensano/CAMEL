import unittest


from camel.app.core.command import Command
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.checkv.checkv import CheckV
from camel.tests import resourceIntensiveTest, longRunningTest


class TestCheckV(CamelTestSuite):
    """
    Tests the CheckV tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('checkv')
    input_fasta = test_file_dir / 'contigs_hev.fasta'

    def test_dependencies(self) -> None:
        """
        Tests if the tool dependencies are available.
        :return: None
        """
        checkm = CheckV()
        for dependency in checkm.dependencies:
            command = Command(f'module load {dependency};')
            command.run(self.running_dir)
            self.assertEqual(command.returncode, 0, f"Dependency '{dependency}' cannot be loaded")

    @longRunningTest()
    @resourceIntensiveTest(reason='RAM usage')
    def test_checkv_tool(self) -> None:
        """
        Tests the CheckV tool.
        :return: None
        """
        checkv = CheckV()
        checkv.add_input_files({'FASTA': [ToolIOFile(TestCheckV.input_fasta)]})
        checkv.update_parameters(threads=4)
        checkv.run(self.running_dir)
        self.assertIn('TSV_quality_summary', checkv.tool_outputs)


if __name__ == '__main__':
    unittest.main()
