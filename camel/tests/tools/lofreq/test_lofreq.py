import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.components.vcf.vcfutils import VCFUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.lofreq.lofreqcall import LofreqCall


class TestLofreq(CamelTestSuite):
    """
    Initializes this testing tool
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('lofreq')
    FILE_FASTA = ToolIOFile(test_file_dir / 'H37Rv.fasta')
    FILE_BAM = ToolIOFile(test_file_dir / 'aln_subsampled.bam')

    def test_lofreq_call(self) -> None:
        """
        Testing Lofreq call on illumina sequencing data
        :return: None
        """
        lofreq = LofreqCall()
        lofreq.add_input_files({'FASTA': [TestLofreq.FILE_FASTA], 'BAM': [TestLofreq.FILE_BAM]})
        lofreq.run(self.running_dir)
        self.verify_output_files(lofreq, 'VCF')
        self.assertGreater(VCFUtils.count_variants(lofreq.tool_outputs['VCF'][0].path), 0)


if __name__ == '__main__':
    unittest.main()
