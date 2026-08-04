import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.art.art import ART


class TestART(CamelTestSuite):
    """
    Tests the ART tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('art')
    input_fasta = test_file_dir / 'example_contigs.fasta'

    def test_art(self) -> None:
        """
        Tests the ART tool.
        :return: None
        """
        art = ART()
        art.add_input_files({
            'FASTA': [ToolIOFile(Path(TestART.input_fasta))]
        })
        art.run(self.running_dir)
        self.verify_output_files(art, 'FASTQ_PE', 2)

    def test_art_adapted_parameters(self) -> None:
        """
        Tests the ART tool with adapted parameters.
        :return: None
        """
        art = ART()
        art.add_input_files({
            'FASTA': [ToolIOFile(Path(TestART.input_fasta))]
        })
        art.update_parameters(read_length=100, fold_coverage=15, mean_size=150, standard_deviation=5)
        art.run(self.running_dir)
        self.verify_output_files(art, 'FASTQ_PE', 2)


if __name__ == '__main__':
    unittest.main()
