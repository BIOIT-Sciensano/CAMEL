import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.gatk4.gatk4applybqsr import GATK4ApplyBQSR
from camel.app.tools.gatk4.gatk4applyvqsr import GATK4ApplyVQSR
from camel.app.tools.gatk4.gatk4baserecalibrator import GATK4BaseRecalibrator
from camel.app.tools.gatk4.gatk4combinegvcfs import GATK4CombineGVCFs
from camel.app.tools.gatk4.gatk4fastaalternatereferencemaker import (
    GATK4FastaAlternateReferenceMaker,
)
from camel.app.tools.gatk4.gatk4gatherbqsrreports import GATK4GatherBQSRReports
from camel.app.tools.gatk4.gatk4genotypegvcfs import GATK4GenotypeGVCFs
from camel.app.tools.gatk4.gatk4haplotypecaller import GATK4HaplotypeCaller
from camel.app.tools.gatk4.gatk4indexfeaturefile import GATK4IndexFeatureFile
from camel.app.tools.gatk4.gatk4selectvariants import GATK4SelectVariants
from camel.app.tools.gatk4.gatk4validatevariants import GATK4ValidateVariants
from camel.app.tools.gatk4.gatk4variantfiltration import GATK4VariantFiltration
from camel.app.tools.gatk4.gatk4variantrecalibrator import GATK4VariantRecalibrator


class TestGATK4(CamelTestSuite):
    """
    Tests the GATK4 tool suite.
    """
    # Get test file and reference file directories
    test_file_dir = CamelTestSuite.get_test_file_dir('gatk4')
    ref_file_dir = CamelTestSuite.get_reference_file_dir('Human', 'GATK-BroadIns', 'hg38', 'v0')

    # Create ToolIOFile input files
    FILE_BAM = ToolIOFile(test_file_dir / 'aln_rg.bam')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'reference.fasta')
    FILE_BQSR = ToolIOFile(test_file_dir / "recal_data.table")
    FILE_KNOWN_SITES = ToolIOFile(test_file_dir / "known_sites.vcf.gz")
    FILE_TXT_INTERVALS = ToolIOFile(test_file_dir / "interval_file.intervals")
    FILE_VCF_JOINTGT = ToolIOFile(test_file_dir / "joint_gt_chr22.vcf.gz")
    FILE_FASTA_REF_hg = ToolIOFile(ref_file_dir / "Homo_sapiens_assembly38.fasta")
    FILE_TXT_RECAL1 = ToolIOFile(test_file_dir / "1_recal_data.csv")
    FILE_TXT_RECAL2 = ToolIOFile(test_file_dir / "2_recal_data.csv")

    # initialize variables
    FILE_gVCF1 = test_file_dir / "NA12877_sub1.g.vcf.gz"
    FILE_gVCF2 = test_file_dir / "NA12877_sub2.g.vcf.gz"
    FILE_VCF = test_file_dir / "var1.vcf"

    def test_gatk4_applybqsr(self) -> None:
        """
        Test GATK4ApplyBQSR
        :return: None
        """
        apply_bqsr = GATK4ApplyBQSR()
        apply_bqsr.add_input_files({
            'BAM': [TestGATK4.FILE_BAM],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF],
            'BQSR': [TestGATK4.FILE_BQSR]
        })
        apply_bqsr.run(self.running_dir)
        self.verify_output_files(apply_bqsr, 'BAM')

    def test_gatk4_applyvqsr(self) -> None:
        """
        Test GATK4ApplyVQSR
        :return: None
        """
        apply_vqsr = GATK4ApplyVQSR()
        apply_vqsr.add_input_files({
            'VCF': [ToolIOFile(TestGATK4.test_file_dir / "joint_gt_chr22.vcf.gz")],
            'TXT_RecalibrationTable': [ToolIOFile(TestGATK4.test_file_dir / "variant_recalibration.tabl")],
            'TXT_tranches': [ToolIOFile(TestGATK4.test_file_dir / "variant_recalibration.tranches")]
        })
        apply_vqsr.update_parameters(
            mode="BOTH",
            filter_level=99.9
        )
        apply_vqsr.run(self.running_dir)
        self.verify_output_files(apply_vqsr, 'VCF')

    def test_gatk4_applyvqsr_err(self) -> None:
        """
        Test GATK4ApplyVQSR tranches error
        :return: None
        """
        apply_vqsr = GATK4ApplyVQSR()
        apply_vqsr.add_input_files({
            'VCF': [ToolIOFile(TestGATK4.test_file_dir / "joint_gt_chr22.vcf.gz")],
            'TXT_RecalibrationTable': [ToolIOFile(TestGATK4.test_file_dir / "variant_recalibration.tabl")],
        })
        apply_vqsr.update_parameters(
            mode="BOTH",
            filter_level=99.9
        )

        with self.assertRaises(InvalidToolInputError):
            apply_vqsr.run(self.running_dir)

    def test_gatk4_baserecalibrator(self) -> None:
        """
        Test GATK4BaseRecalibrator
        :return: None
        """
        baserecalibrator = GATK4BaseRecalibrator()
        baserecalibrator.add_input_files({
            'BAM': [TestGATK4.FILE_BAM],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF],
            'VCF_KNOWN_SNPS': [TestGATK4.FILE_KNOWN_SITES]
        })
        baserecalibrator.run(self.running_dir)
        self.verify_output_files(baserecalibrator, 'TXT_RecalibrationTable')

    def test_gatk4_combinegvcf(self) -> None:
        """
        Test GATK4CombineGVCFs
        :return: None
        """
        combinegvcf = GATK4CombineGVCFs()
        combinegvcf.add_input_files({
            'gVCF': [ToolIOFile(TestGATK4.FILE_gVCF1), ToolIOFile(TestGATK4.FILE_gVCF2)],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF_hg],
        })
        combinegvcf.run(self.running_dir)
        self.verify_output_files(combinegvcf, 'gVCF')

    def test_gatk4_fastaalternatereferencemaker(self) -> None:
        """
        Test GATK4FastaAlternateReferenceMaker
        :return: None
        """
        fastaalternatereferencemaker = GATK4FastaAlternateReferenceMaker()
        fastaalternatereferencemaker.add_input_files({
            'VCF': [TestGATK4.FILE_KNOWN_SITES],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF],
        })
        fastaalternatereferencemaker.update_parameters(**{'concatenate_sequence_segments': 'true'})
        fastaalternatereferencemaker.run(self.running_dir)
        self.verify_output_files(fastaalternatereferencemaker, 'FASTA')

    def test_gatk4_gatherbqsrreports(self) -> None:
        """
        Test GATK4GatherBQSRReports
        :return: None
        """
        gatherbqsrreports = GATK4GatherBQSRReports()
        gatherbqsrreports.add_input_files({
            'TXT_intervals': [TestGATK4.FILE_TXT_RECAL1, TestGATK4.FILE_TXT_RECAL2],
        })
        gatherbqsrreports.run(self.running_dir)
        self.verify_output_files(gatherbqsrreports, 'TXT_RecalibrationTable')

    def test_gatk4_genotypegvcfs(self) -> None:
        """
        Test GATK4GenotypeGVCFs
        :return: None
        """
        genotypegvcfs = GATK4GenotypeGVCFs()
        genotypegvcfs.add_input_files({
            'gVCF': [ToolIOFile(TestGATK4.FILE_gVCF1)],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF_hg],
        })
        genotypegvcfs.run(self.running_dir)
        self.verify_output_files(genotypegvcfs, 'VCF_MultipleSample')

    def test_gatk4_haplotypecaller(self) -> None:
        """
        Test GATK4HaplotypeCaller
        :return: None
        """
        haplotypecaller = GATK4HaplotypeCaller()
        haplotypecaller.add_input_files({
            'BAM': [TestGATK4.FILE_BAM],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF],
        })
        haplotypecaller.run(self.running_dir)
        self.verify_output_files(haplotypecaller, 'VCF')

    def test_gatk4_indexfeaturefile(self) -> None:
        """
        Test GATK4IndexFeatureFile
        :return: None
        """
        indexfeaturefile = GATK4IndexFeatureFile()

        # Output automatically made in directory of input file
        # Symlink to input file in running_dir
        vcf_workingdir = Path(self.running_dir) / "test_indexFeatureFile.vcf"
        vcf_original = TestGATK4.FILE_VCF
        vcf_workingdir.symlink_to(vcf_original)

        indexfeaturefile.add_input_files({
            'VCF': [ToolIOFile(vcf_workingdir)]
        })
        indexfeaturefile.run(self.running_dir)
        self.verify_output_files(indexfeaturefile, 'IDX')

    def test_gatk4_indexfeaturefile_gz(self) -> None:
        """
        Tests the GATK4IndexFeatureFile with GZIP.
        :return: None
        """
        indexfeaturefile_gz = GATK4IndexFeatureFile()

        # Output automatically made in directory of input file
        # Symlink to input file in running_dir
        vcf_workingdir = Path(self.running_dir) / "test_indexFeatureFile.vcf.gz"
        vcf_original = TestGATK4.FILE_gVCF1
        vcf_workingdir.symlink_to(vcf_original)

        indexfeaturefile_gz.add_input_files({
            'VCF_gz': [ToolIOFile(vcf_workingdir)]
        })
        indexfeaturefile_gz.run(self.running_dir)
        self.verify_output_files(indexfeaturefile_gz, 'IDX')

    def test_gatk4_selectvariants(self) -> None:
        """
        Test GATK4SelectVariants
        :return: None
        """
        selectvariants = GATK4SelectVariants()
        selectvariants.add_input_files({
            'VCF': [ToolIOFile(TestGATK4.FILE_VCF)],
        })
        selectvariants.run(self.running_dir)
        self.verify_output_files(selectvariants, 'VCF')

    def test_gatk4_validatevariants(self) -> None:
        """
        Test GATK4ValidateVariants
        :return: None
        """
        validatevariants = GATK4ValidateVariants()
        validatevariants.add_input_files({
            'VCF': [ToolIOFile(TestGATK4.FILE_VCF)],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF],
        })
        validatevariants.run(self.running_dir)
        self.verify_output_files(validatevariants, 'TXT_metrics')

    def test_gatk4_variantfiltration(self) -> None:
        """
        Test GATK4VariantFiltration
        :return: None
        """
        variantfiltration = GATK4VariantFiltration()
        variantfiltration.add_input_files({
            'VCF': [ToolIOFile(TestGATK4.FILE_VCF)],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF],
        })
        variantfiltration.update_parameters(**{'filter-names': 'HardQC_filter,lowDP',
                                'filter-expressions': 'QD<2.0||FS>60.0||MQ<40.0||ReadPosRankSum<-8.0,DP<200'})
        variantfiltration.run(self.running_dir)
        self.verify_output_files(variantfiltration, 'VCF')

    def test_gatk4_variantrecalibrator(self) -> None:
        """
        Test GATK4VariantRecalibrator
        :return: None
        """
        variantrecalibrator = GATK4VariantRecalibrator()
        variantrecalibrator.add_input_files({
            'VCF': [TestGATK4.FILE_VCF_JOINTGT],
            'FASTA_REF': [TestGATK4.FILE_FASTA_REF_hg]
        })
        variantrecalibrator.update_parameters(
            resources = f"hapmap,known=false,training=true,truth=true,prior=15.0,{TestGATK4.ref_file_dir}/hapmap_3.3.hg38.vcf.gz",
            use_annotation = "DP",
            mode = "BOTH"
        )
        variantrecalibrator.run(self.running_dir)
        self.verify_output_files(variantrecalibrator, 'TXT_RecalibrationTable')
        self.verify_output_files(variantrecalibrator, 'TXT_tranches')

if __name__ == '__main__':
    unittest.main()
