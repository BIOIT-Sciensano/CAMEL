import unittest

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.quast.quast import Quast


class TestQuast(CamelTestSuite):
    """
    Tests the Quast tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('quast')
    FASTA_IN = ToolIOFile(test_file_dir / 'contigs_listeria.fasta')
    FASTA_REF = ToolIOFile(test_file_dir / 'ref_listeria.fasta')
    BAM_IN = ToolIOFile(test_file_dir / 'output_contigs.bam')
    BAM_REF = ToolIOFile(test_file_dir / 'output_ref.bam')

    def test_quast(self) -> None:
        """
        Tests Quast with default options.
        :return: None
        """
        quast = Quast()
        quast.add_input_files({'FASTA': [TestQuast.FASTA_IN]})
        quast.run(self.running_dir)
        self.verify_output_files(quast, 'TSV')
        self.verify_output_files(quast, 'HTML')

    def test_quast_with_ref(self) -> None:
        """
        Tests Quast with a reference genome.
        :return: None
        """
        quast = Quast()
        quast.add_input_files({'FASTA': [TestQuast.FASTA_IN], 'FASTA_Ref': [TestQuast.FASTA_REF]})
        quast.run(self.running_dir)
        self.verify_output_files(quast, 'TSV')
        self.verify_output_files(quast, 'HTML')

    def test_quast_with_glimmer(self) -> None:
        """
        Tests Quast with glimmer gene prediction enabled.
        :return: None
        """
        quast = Quast()
        quast.add_input_files({'FASTA': [TestQuast.FASTA_IN], 'FASTA_Ref': [TestQuast.FASTA_REF]})
        quast.update_parameters(glimmer=True)
        quast.run(self.running_dir)
        self.verify_output_files(quast, 'TSV')
        self.verify_output_files(quast, 'HTML')
        self.verify_output_files(quast, 'GFF')

    def test_quast_with_bam_input(self) -> None:
        """
        Tests Quast with a BAM input file containing reads mapped against contigs.
        :return: None
        """
        quast = Quast()
        quast.add_input_files({'FASTA': [TestQuast.FASTA_IN], 'BAM': [TestQuast.BAM_IN]})
        quast.run(self.running_dir)
        self.verify_output_files(quast, 'TSV')
        self.verify_output_files(quast, 'HTML')

    def test_quast_with_bam_input_and_ref(self) -> None:
        """
        Tests Quast with a BAM input file containing reads mapped against contigs.
        :return: None
        """
        quast = Quast()
        quast.add_input_files({
            'FASTA': [TestQuast.FASTA_IN],
            'FASTA_Ref': [TestQuast.FASTA_REF],
            'BAM': [TestQuast.BAM_IN],
            'BAM_Ref': [TestQuast.BAM_REF]
        })
        quast.run(self.running_dir)
        self.verify_output_files(quast, 'TSV')
        self.verify_output_files(quast, 'HTML')


if __name__ == '__main__':
    unittest.main()
