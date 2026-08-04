import unittest

from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import vcfutils

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.lofreq.lofreqcall import LofreqCall
from camel.app.tools.lofreq.lofreqindelqual import LofreqIndelqual


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
        self.assertGreater(vcfutils.count_variants(lofreq.tool_outputs['VCF'][0].path), 0)

    def test_lofreq_indelqual(self) -> None:
        """
        Testing Lofreq indelqual command on illumina BAM file
        :return: None
        """
        lofreq_indelqual = LofreqIndelqual()
        lofreq_indelqual.add_input_files({'FASTA': [TestLofreq.FILE_FASTA], 'BAM': [TestLofreq.FILE_BAM]})
        lofreq_indelqual.run(self.running_dir)
        self.verify_output_files(lofreq_indelqual, 'BAM')


if __name__ == '__main__':
    unittest.main()
