import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.mash.mashscreen import MashScreen


class TestMash(CamelTestSuite):
    """
    Tests the mash tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mash')
    input_db = test_file_dir / 'influenza_a-subtypes.msh'
    input_fastq = [test_file_dir / 'influenza_a_1.fastq.gz', test_file_dir / 'influenza_a_2.fastq.gz']
    input_fasta = test_file_dir / 'influenza_a-full_genome.fasta'

    def test_mash_screen(self) -> None:
        """
        Tests the mash screen tool.
        :return: None
        """
        mash_screen = MashScreen()
        mash_screen.add_input_files({
            'FASTQ': [ToolIOFile(x) for x in TestMash.input_fastq],
            'DB': [ToolIOFile(TestMash.input_db)]
        })
        mash_screen.run(self.running_dir)
        self.verify_output_files(mash_screen, 'TSV')

    def test_mash_screen_fasta_input(self) -> None:
        """
        Tests the mash screen tool with FASTA input.
        :return: None
        """
        mash_screen = MashScreen()
        mash_screen.add_input_files({
            'FASTA': [ToolIOFile(TestMash.input_fasta)],
            'DB': [ToolIOFile(TestMash.input_db)]
        })
        mash_screen.run(self.running_dir)
        self.verify_output_files(mash_screen, 'TSV')

    def test_mash_screen_params(self) -> None:
        """
        Tests the mash screen tool with updated parameters.
        :return: None
        """
        mash_screen = MashScreen()
        mash_screen.add_input_files({
            'FASTQ': [ToolIOFile(x) for x in TestMash.input_fastq],
            'DB': [ToolIOFile(TestMash.input_db)]
        })
        mash_screen.update_parameters(max_p_value=0.05, min_identity=0.5, winner_takes_all=True, threads=2)
        mash_screen.run(self.running_dir)
        self.verify_output_files(mash_screen, 'TSV')


if __name__ == '__main__':
    unittest.main()
