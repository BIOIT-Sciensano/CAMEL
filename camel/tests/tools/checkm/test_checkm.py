import unittest
from pathlib import Path


from camel.app.core.command import Command
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.checkm.checkm import CheckM
from camel.tests import longRunningTest, resourceIntensiveTest


class TestCheckM(CamelTestSuite):
    """
    Tests the CheckM tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('checkm')
    input_fasta = test_file_dir / 'contigs_neisseria.fasta'

    def test_dependencies(self) -> None:
        """
        Tests if the tool dependencies are available.
        :return: None
        """
        checkm = CheckM()
        for dependency in checkm.dependencies:
            command = Command(f'module load {dependency};')
            command.run(self.running_dir)
            self.assertEqual(command.returncode, 0, f"Dependency '{dependency}' cannot be loaded")

    @longRunningTest()
    @resourceIntensiveTest(reason='RAM usage')
    def test_checkm(self) -> None:
        """
        Tests the CheckM tool.
        :return: None
        """
        checkm = CheckM()
        checkm.add_input_files({'FASTA': [ToolIOFile(TestCheckM.input_fasta)]})
        checkm.update_parameters(reduced_tree=True)
        checkm.run(self.running_dir)
        self.assertIn('TSV', checkm.tool_outputs)
        self.assertGreater(Path(checkm.tool_outputs['TSV'][0].path).stat().st_size, 0)
        self.assertIn('results', checkm.informs)


if __name__ == '__main__':
    unittest.main()
