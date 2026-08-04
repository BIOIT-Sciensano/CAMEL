import unittest

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.scripts.variantcalling.medaka.maincallingmedaka import MainCallingMedaka


class TestMedaka(CamelTestSuite):
    """
    Tests the Medaka tool suite.
    """
    # Get test files and reference files directories
    test_file_dir = CamelTestSuite.get_test_file_dir('medaka')

    # Create ToolIOFile input files
    FILE_BAM = ToolIOFile(test_file_dir / 'calls_to_draft_subsampled.bam')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'contig_1.fasta')

    def test_medaka_main_variant_calling(self) -> None:
        """
        Tests the variant calling main script.
        :return: None
        """
        output_file_vcf = self.running_dir / 'calls_to_draft_subsampled.vcf'
        args = [
            '--bam',
            str(TestMedaka.FILE_BAM),
            '--reference',
            str(TestMedaka.FILE_FASTA_REF),
            '--working-dir',
            str(self.running_dir),
            '--output',
            str(output_file_vcf),
        ]
        main_calling_variant = MainCallingMedaka(args)
        main_calling_variant.run()
        self.assertGreater(output_file_vcf.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
