import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.utils import vcfutils
from camel.app.tools.sniffles.sniffles import Sniffles


class TestSniffles(CamelTestSuite):
    """
    Initializes the Sniffles testing tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('clair3')
    FILE_FASTA = ToolIOFile(test_file_dir / 'bsubtilis.fa')
    FILE_BAM_ONT = ToolIOFile(test_file_dir / 'bsubtilis_ont.bam')

    def test_sniffles(self) -> None:
        """
        Actually testing Sniffles on ONT sequencing data.
        """
        sniffles = Sniffles()
        sniffles.add_input_files({'BAM': [TestSniffles.FILE_BAM_ONT], 'FASTA': [TestSniffles.FILE_FASTA]})
        sniffles.run(self.running_dir)
        self.verify_output_files(sniffles, 'VCF')
        self.assertGreater(vcfutils.count_variants(sniffles.tool_outputs['VCF'][0].path), 0)


if __name__ == '__main__':
    unittest.main()
