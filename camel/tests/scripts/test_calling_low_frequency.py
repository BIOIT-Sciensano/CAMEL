import traceback
import unittest

from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.utils import vcfutils

from camel.scripts.variantcalling.lofreq.maincallinglofreq import main


class TestVariantCalling(CamelTestSuite):
    """
    Tests the variant calling tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('variant_calling_lofreq')
    input_fasta_ref_file = test_file_dir / 'ref.fasta'
    input_gff_file = test_file_dir / 'ref.gff'
    fastq_1 = test_file_dir / 'read_R1.fq.gz'
    fastq_2 = test_file_dir / 'read_R2.fq.gz'

    def test_variant_calling(self) -> None:
        """
        Tests the variant calling main script.
        :return: None
        """
        output_file_vcf = self.running_dir / 'detected_variants.vcf'
        result = cliutils.invoke(
            main, [
                '--fastq-pe', str(TestVariantCalling.fastq_1), str(TestVariantCalling.fastq_2),
                '--reference', str(TestVariantCalling.input_fasta_ref_file),
                '--gff', str(TestVariantCalling.input_gff_file),
                '--working-dir', str(self.running_dir),
                '--output-vcf', str(output_file_vcf),
                '--input-type', 'illumina',
                '--call-indels'
            ]
        )
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_file_vcf.stat().st_size, 0)
        self.assertGreater(vcfutils.count_variants(output_file_vcf), 0)
