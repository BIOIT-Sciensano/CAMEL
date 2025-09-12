import shutil
import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.error import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.bcftools.bcftoolscall import BcftoolsCall
from camel.app.tools.bcftools.bcftoolsconsensus import BcftoolsConsensus
from camel.app.tools.bcftools.bcftoolscsq import BcftoolsCsq
from camel.app.tools.bcftools.bcftoolsfilter import BcftoolsFilter
from camel.app.tools.bcftools.bcftoolsindex import BcftoolsIndex
from camel.app.tools.bcftools.bcftoolsmpileup import BcftoolsMpileup
from camel.app.tools.bcftools.bcftoolsnorm import BcftoolsNorm
from camel.app.tools.bcftools.bcftoolsview import BcftoolsView


class TestBcftools(CamelTestSuite):
    """
    Tests the bcftools tool suite.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('bcftools')
    FILE_VCF = ToolIOFile(test_file_dir / 'variants.vcf')
    FILE_VCF_GZ = ToolIOFile(test_file_dir / 'variants.vcf.gz')
    FILE_PILEUP = ToolIOFile(test_file_dir / 'pileup.vcf')
    FILE_PILEUP_GZ = ToolIOFile(test_file_dir / 'pileup.vcf.gz')
    FILE_BED = ToolIOFile(test_file_dir / 'regions.bed')
    FILE_FASTA = ToolIOFile(test_file_dir / 'reference_h37Rv.fasta')
    FILE_GFF = ToolIOFile(test_file_dir / 'annotation_h37Rv.gff')
    FILE_BAM_TOY = ToolIOFile(test_file_dir / 'toy.bam')
    FILE_FASTA_TOY = ToolIOFile(test_file_dir / 'toy.fasta')

    def test_bcftools_call(self) -> None:
        """
        Tests the bcftools call tool.
        :return: None
        """
        bcftools_call = BcftoolsCall()
        bcftools_call.add_input_files({'VCF': [TestBcftools.FILE_PILEUP]})
        bcftools_call.update_parameters(ploidy=1, output_type='v')
        bcftools_call.run(self.running_dir)
        self.verify_output_files(bcftools_call, 'VCF')

    def test_bcftools_call_gz_input(self) -> None:
        """
        Tests the bcftools call tool with gzipped input.
        :return: None
        """
        bcftools_call = BcftoolsCall()
        bcftools_call.add_input_files({'VCF_GZ': [TestBcftools.FILE_PILEUP_GZ]})
        bcftools_call.update_parameters(ploidy=1, output_type='v')
        bcftools_call.run(self.running_dir)
        self.verify_output_files(bcftools_call, 'VCF')

    def test_bcftools_call_gz_output(self) -> None:
        """
        Tests the bcftools call tool with gzipped input.
        :return: None
        """
        bcftools_call = BcftoolsCall()
        bcftools_call.add_input_files({'VCF_GZ': [TestBcftools.FILE_PILEUP_GZ]})
        bcftools_call.update_parameters(ploidy=1, output_type='z', output_filename='calls.vcf.gz')
        bcftools_call.run(self.running_dir)
        self.verify_output_files(bcftools_call, 'VCF_GZ')

    def test_bcftools_consensus(self) -> None:
        """
        Tests the bcftools consensus tool.
        :return: None
        """
        bcftools_consensus = BcftoolsConsensus()
        bcftools_consensus.add_input_files({
            'VCF_GZ': [TestBcftools.FILE_VCF_GZ],
            'FASTA': [TestBcftools.FILE_FASTA]
        })
        bcftools_consensus.run(self.running_dir)
        self.verify_output_files(bcftools_consensus, 'FASTA')

    def test_bcftools_csq(self) -> None:
        """
        Tests bcftools csq.
        :return: None
        """
        bcftools_csq = BcftoolsCsq()
        bcftools_csq.add_input_files({
            'VCF_GZ': [TestBcftools.FILE_VCF_GZ],
            'FASTA': [TestBcftools.FILE_FASTA],
            'GFF': [TestBcftools.FILE_GFF]
        })
        bcftools_csq.run(self.running_dir)
        self.verify_output_files(bcftools_csq, 'VCF')

    def test_bcftools_filter(self) -> None:
        """
        Tests bcftools filter with valid input.
        :return: None
        """
        bcftools_filter = BcftoolsFilter()
        bcftools_filter.add_input_files({
            'VCF_GZ': [TestBcftools.FILE_VCF_GZ],
            'BED': [TestBcftools.FILE_BED]
        })
        bcftools_filter.run(self.running_dir)
        self.verify_output_files(bcftools_filter, 'VCF')

    def test_bcftools_filter_invalid_input(self) -> None:
        """
        Tests that bcftools filter fails with invalid input.
        :return: None
        """
        bcftools_filter = BcftoolsFilter()
        bcftools_filter.add_input_files({
            'VCF_GZ': [TestBcftools.FILE_FASTA],
            'BED': [TestBcftools.FILE_BED]
        })
        with self.assertRaises(ToolExecutionError):
            bcftools_filter.run(self.running_dir)

    def test_bcftools_index(self) -> None:
        """
        Tests the bcftools index function.
        :return: None
        """
        bcftools_index = BcftoolsIndex()
        bcftools_index.add_input_files({'VCF_GZ': [TestBcftools.FILE_VCF_GZ]})
        bcftools_index.run(self.running_dir)
        self.assertTrue(self.running_dir / f"{bcftools_index.tool_outputs['VCF_GZ'][0].path}.csi", 'CSI index not found')

    def test_bcftools_index_no_symlink(self) -> None:
        """
        Tests the bcftools index function without creating a symlink.
        :return: None
        """
        # Copy the input VCF_GZ file to the working directory
        vcf_gz_in = self.running_dir / 'input.vcf.gz'
        shutil.copyfile(str(TestBcftools.FILE_VCF_GZ), str(vcf_gz_in))

        # Create index
        bcftools_index = BcftoolsIndex()
        bcftools_index.add_input_files({'VCF_GZ': [ToolIOFile(vcf_gz_in)]})
        bcftools_index.update_parameters(symlink_input=False)
        bcftools_index.run(self.running_dir)

        # Check if the index was created in the location of the input file
        self.assertTrue(self.running_dir / f'{vcf_gz_in}.csi', 'CSI index not found')

    def test_bcftools_mpileup(self) -> None:
        """
        Tests bcftools mpileup tool.
        :return: None
        """
        bcftools_mpileup = BcftoolsMpileup()
        bcftools_mpileup.add_input_files({'BAM': [TestBcftools.FILE_BAM_TOY], 'FASTA': [TestBcftools.FILE_FASTA_TOY]})
        bcftools_mpileup.run(self.running_dir)
        self.verify_output_files(bcftools_mpileup, 'VCF')

    def test_bcftools_mpileup_gz_out(self) -> None:
        """
        Tests bcftools mpileup tool with VCF_GZ output.
        :return: None
        """
        bcftools_mpileup = BcftoolsMpileup()
        bcftools_mpileup.add_input_files({'BAM': [TestBcftools.FILE_BAM_TOY], 'FASTA': [TestBcftools.FILE_FASTA_TOY]})
        bcftools_mpileup.update_parameters(output_type='z', output_filename='pileup.vcf.gz')
        bcftools_mpileup.run(self.running_dir)
        self.verify_output_files(bcftools_mpileup, 'VCF_GZ')

    def test_bcftools_norm(self) -> None:
        """
        Tests the bcftools norm command.
        :return: None
        """
        bcftools_norm = BcftoolsNorm()
        bcftools_norm.add_input_files({
            'VCF_GZ': [TestBcftools.FILE_VCF_GZ],
            'FASTA': [TestBcftools.FILE_FASTA]
        })
        bcftools_norm.run(self.running_dir)
        self.verify_output_files(bcftools_norm, 'VCF')
        self.assertIn('total', bcftools_norm.informs)

    def test_bcftools_norm_gz_output(self) -> None:
        """
        Tests the bcftools norm command with compressed output.
        :return: None
        """
        bcftools_norm = BcftoolsNorm()
        bcftools_norm.add_input_files({
            'VCF_GZ': [TestBcftools.FILE_VCF_GZ],
            'FASTA': [TestBcftools.FILE_FASTA]
        })
        bcftools_norm.update_parameters(output_type='z', rm_dup='snps')
        bcftools_norm.run(self.running_dir)
        self.verify_output_files(bcftools_norm, 'VCF_GZ')
        self.assertIn('total', bcftools_norm.informs)

    def test_bcftools_view_vcf_to_vcf_gz(self) -> None:
        """
        Tests the bcftools view class to convert a VCF_GZ file to VCF format.
        :return: None
        """
        bcftools_view = BcftoolsView()
        bcftools_view.add_input_files({'VCF': [TestBcftools.FILE_VCF]})
        bcftools_view.update_parameters(output_type='z', output_filename='variants.vcf.gz')
        bcftools_view.run(self.running_dir)
        self.verify_output_files(bcftools_view, 'VCF_GZ')

    def test_bcftools_view_vcf_gz_to_vcf(self) -> None:
        """
        Tests the bcftools view class to convert a VCF_GZ file to VCF format.
        :return: None
        """
        bcftools_view = BcftoolsView()
        bcftools_view.add_input_files({'VCF_GZ': [TestBcftools.FILE_VCF_GZ]})
        bcftools_view.update_parameters(output_type='v')
        bcftools_view.run(self.running_dir)
        self.verify_output_files(bcftools_view, 'VCF')


if __name__ == '__main__':
    unittest.main()
