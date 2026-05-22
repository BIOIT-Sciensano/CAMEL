import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import vcfutils

from camel.app.config import config
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.clair3.clair3 import Clair3


class TestClair3(CamelTestSuite):
    """
    Initializes this testing tool
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('clair3')
    FILE_FASTA = ToolIOFile(test_file_dir / 'bsubtilis.fa')
    FILE_BAM_ILLUMINA = ToolIOFile(test_file_dir / 'bsubtilis_illumina.bam')
    FILE_BAM_ONT = ToolIOFile(test_file_dir / 'bsubtilis_ont.bam')

    def test_clair3_illumina(self) -> None:
        """
        Tests Clair3 on Illumina input data.
        :return: None
        """
        clair3 = Clair3()
        clair3.add_input_files({'FASTA': [TestClair3.FILE_FASTA], 'BAM': [TestClair3.FILE_BAM_ILLUMINA]})
        clair3.update_parameters(output_path=self.running_dir, haploid_precise=True, no_phasing=True, include_ctgs=True)
        clair3.run(self.running_dir)
        self.verify_output_files(clair3, 'VCF')
        self.assertGreater(vcfutils.count_variants(clair3.tool_outputs['VCF'][0].path), 0)

    def test_clair3_ont(self) -> None:
        """
        Tests Clair3 on ONT input data.
        :return: None
        """
        clair3 = Clair3()
        clair3.add_input_files({'FASTA': [TestClair3.FILE_FASTA], 'BAM': [TestClair3.FILE_BAM_ONT]})
        clair3.update_parameters(
            platform='ont',
            model_path=str(Path(config.dir_db, 'clair3', 'models', 'ont')),
            include_ctgs=True
        )

        clair3.run(self.running_dir)
        self.verify_output_files(clair3, 'VCF')
        self.assertGreater(vcfutils.count_variants(clair3.tool_outputs['VCF'][0].path), 0)


if __name__ == '__main__':
    unittest.main()
