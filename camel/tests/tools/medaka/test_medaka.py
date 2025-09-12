import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.medaka.medakainference import MedakaInference
from camel.app.tools.medaka.medakasequence import MedakaSequence
from camel.app.tools.medaka.medakavcf import MedakaVcf


class TestMedaka(CamelTestSuite):
    """
    Tests the Medaka tool suite.
    """
    # Get test files and reference files directories
    test_file_dir = CamelTestSuite.get_test_file_dir('medaka')

    # Create ToolIOFile input files
    FILE_BAM = ToolIOFile(test_file_dir / 'calls_to_draft_subsampled.bam')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'contig_1.fasta')
    FILE_HDF = ToolIOFile(test_file_dir / 'consensus_probs_subsample.hdf')

    def test_medaka_inference(self) -> None:
        """
        Tests MedakaInference.
        :return: None
        """
        inference = MedakaInference()
        inference.add_input_files({'BAM': [TestMedaka.FILE_BAM]})
        inference.run(self.running_dir)
        self.verify_output_files(inference, 'HDF')

    def test_medaka_sequence(self) -> None:
        """
        Tests MedakaSequence.
        :return: None
        """
        sequence = MedakaSequence()
        sequence.add_input_files({
            'HDF': [TestMedaka.FILE_HDF],
            'FASTA': [TestMedaka.FILE_FASTA_REF]
        })
        sequence.run(self.running_dir)
        self.verify_output_files(sequence, 'FASTA')

    def test_medaka_vcf(self) -> None:
        """
        Tests MedakaVcf.
        :return: None
        """
        vcf = MedakaVcf()
        vcf.add_input_files({
            'HDF': [TestMedaka.FILE_HDF],
            'FASTA': [TestMedaka.FILE_FASTA_REF]
        })
        vcf.run(self.running_dir)
        self.verify_output_files(vcf, 'VCF')


if __name__ == '__main__':
    unittest.main()
