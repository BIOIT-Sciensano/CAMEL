import unittest


from camel.app.core.command import Command
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.busco.busco import Busco


class TestBusco(CamelTestSuite):
    """
    Tests the BUSCO tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('checkm')
    input_fasta = test_file_dir / 'contigs_neisseria.fasta'
    input_fasta_2 = test_file_dir / 'busco_failure_noble.fasta'

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

    def test_busco_failure_noble(self) -> None:
        """
        Tests the busco tool on a fasta file that initially failed in noble and was not caught by the previous
        'test_busco' unittest. The issue was resolved on noble by replacing the jammy tarball with the prebuilt binary
        (https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux).
        :return: None
        """
        busco = Busco()
        busco.add_input_files({'FASTA': [ToolIOFile(TestBusco.input_fasta_2)]})
        busco.update_parameters(lineage_dataset='bacteria_odb10')
        busco.run(self.running_dir)
        self.verify_output_files(busco, 'TXT')
        self.assertIn('results', busco.informs)


if __name__ == '__main__':
    unittest.main()
