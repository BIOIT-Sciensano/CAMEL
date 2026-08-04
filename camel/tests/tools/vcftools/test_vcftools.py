import unittest

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.vcftools.vcftoolsannotate import VCFtoolsAnnotate


class TestVCFtools(CamelTestSuite):
    """
    Tests the vcftools tool suite.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('bcftools')
    FILE_VCF_GZ = ToolIOFile(test_file_dir / 'variants.vcf.gz')

    def test_vcftools_vcfannotate(self) -> None:
        """
        Tests vcftools vcf-annotate.
        :return: None
        """
        vcf_annotate = VCFtoolsAnnotate()
        vcf_annotate.add_input_files({
            'VCF_GZ': [TestVCFtools.FILE_VCF_GZ]
        })
        vcf_annotate.update_parameters(
            fill_hwe='',
            fill_icf=''
        )
        vcf_annotate.run(self.running_dir)
        self.verify_output_files(vcf_annotate, 'VCF_GZ')


if __name__ == '__main__':
    unittest.main()
