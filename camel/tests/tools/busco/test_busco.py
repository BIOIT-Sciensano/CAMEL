import unittest


from camel.app.command.command import Command
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.busco.busco import Busco


class TestBusco(CamelTestSuite):
    """
    Tests the BUSCO tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('checkm')
    input_fasta = test_file_dir / 'contigs_neisseria.fasta'

    def test_dependencies(self) -> None:
        """
        Tests if the tool dependencies are available.
        :return: None
        """
        busco = Busco()
        for dependency in busco.dependencies:
            command = Command(f'module load {dependency};')
            command.run(self.running_dir)
            self.assertEqual(command.returncode, 0, f"Dependency '{dependency}' cannot be loaded")

    def test_busco(self) -> None:
        """
        Tests the busco tool.
        :return: None
        """
        busco = Busco()
        busco.add_input_files({'FASTA': [ToolIOFile(TestBusco.input_fasta)]})
        busco.update_parameters(lineage_dataset='bacteria_odb10')
        busco.run(self.running_dir)
        self.verify_output_files(busco, 'TXT')
        self.assertIn('results', busco.informs)


if __name__ == '__main__':
    unittest.main()
