import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.bakta.bakta import Bakta
from camel.tests import longRunningTest


class TestBakta(CamelTestSuite):
    """
    Tests the BAKTA tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('bakta')
    input_fasta = test_file_dir / 'E_coli_K12_subset.fasta'
    input_spp = test_file_dir / 'B_subtilis_168_subset.fasta'
    input_meta_fasta = test_file_dir / 'metagenome_subset.fasta'

    @longRunningTest()
    def test_bakta(self) -> None:
        """
        Tests the Bakta tool.
        :return: None
        """
        bakta = Bakta()
        bakta.add_input_files({
            'FASTA': [ToolIOFile(Path(TestBakta.input_fasta))]
        })
        bakta.run(self.running_dir)
        self.verify_output_files(bakta, 'FAA', 1)
        self.verify_output_files(bakta, 'GFF3', 1)
        self.verify_output_files(bakta, 'GBFF', 1)

    @longRunningTest()
    def test_bakta_args(self) -> None:
        """
        Tests the Bakta tool with added arguments.
        :return: None
        """
        bakta = Bakta()
        bakta.add_input_files({
            'FASTA': [ToolIOFile(Path(TestBakta.input_spp))]
        })
        bakta.update_parameters(genus='bacillus', species='subtilis', strain='168')
        bakta.run(self.running_dir)
        self.verify_output_files(bakta, 'FAA', 1)
        self.verify_output_files(bakta, 'GFF3', 1)
        self.verify_output_files(bakta, 'GBFF', 1)

    @longRunningTest()
    def test_bakta_meta(self) -> None:
        """
        Tests the Bakta tool in metagenome mode.
        :return: None
        """
        bakta = Bakta()
        bakta.add_input_files({
            'FASTA': [ToolIOFile(Path(TestBakta.input_meta_fasta))]
        })
        bakta.update_parameters(meta=None)
        bakta.run(self.running_dir)
        self.verify_output_files(bakta, 'FAA', 1)
        self.verify_output_files(bakta, 'GFF3', 1)
        self.verify_output_files(bakta, 'GBFF', 1)


if __name__ == '__main__':
    unittest.main()
