import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import fastautils

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.errors import InvalidParameterError, InvalidToolInputError
from camel.app.tools.gatk.gatkanalyzecovariates import GATKAnalyzeCovariates
from camel.app.tools.gatk.gatkbaserecalibrator import GATKBaseRecalibrator
from camel.app.tools.gatk.gatkfastaalternatereferencemaker import GATKFastaAlternateReferenceMaker
from camel.app.tools.gatk.gatkgenotypegvcfs import GATKGenotypeGVCFs
from camel.app.tools.gatk.gatkhaplotypecaller import GATKHaplotypeCaller
from camel.app.tools.gatk.gatkindelrealigner import GATKIndelRealigner
from camel.app.tools.gatk.gatkmutect2 import GATKMuTect2
from camel.app.tools.gatk.gatkprintreads import GATKPrintReads
from camel.app.tools.gatk.gatkrealignertargetcreator import GATKRealignerTargetCreator
from camel.app.tools.gatk.gatkselectvariants import GATKSelectVariants


class TestGATK(CamelTestSuite):
    """
    Tests the GATK (v3) tools.
    Test data: phiX reads and reference from the Galaxy GATK2 wrappers (galaxyproject/tools-iuc), the variants were
    called with GATK 3.7 HaplotypeCaller (single SNP: phiX174:1424 C>A).
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('gatk')
    FILE_BAM = ToolIOFile(test_file_dir / 'fake_phiX_reads_1.bam')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'phiX.fasta')
    FILE_VCF = ToolIOFile(test_file_dir / 'phiX_variants.vcf')

    def __create_dir(self, name: str) -> Path:
        """
        Creates a subdirectory of the running directory (used to chain tools).
        :param name: Directory name
        :return: Path to the directory
        """
        dir_ = self.running_dir / name
        dir_.mkdir()
        return dir_

    def __run_base_recalibrator(self, bam: ToolIOFile, name: str) -> ToolIOFile:
        """
        Runs the base recalibrator.
        :param bam: Alignment file
        :param name: Name of the running directory
        :return: Recalibration table
        """
        base_recalibrator = GATKBaseRecalibrator()
        base_recalibrator.add_input_files({
            'BAM': [bam],
            'FASTA_REF': [TestGATK.FILE_FASTA_REF],
            'VCF_KNOWN_SNPS': [TestGATK.FILE_VCF]
        })
        base_recalibrator.run(self.__create_dir(name))
        self.verify_output_files(base_recalibrator, 'TXT_RecalibrationTable')
        return base_recalibrator.tool_outputs['TXT_RecalibrationTable'][0]

    def test_version(self) -> None:
        """
        Tests the version retrieval.
        :return: None
        """
        self.assertRegex(GATKHaplotypeCaller().version, r'^3\.\d+')

    def test_realigner_target_creator_indel_realigner(self) -> None:
        """
        Tests the realigner target creator, followed by the indel realigner.
        :return: None
        """
        target_creator = GATKRealignerTargetCreator()
        target_creator.add_input_files({'BAM': [TestGATK.FILE_BAM], 'FASTA_REF': [TestGATK.FILE_FASTA_REF]})
        target_creator.run(self.__create_dir('target_creator'))
        self.assertIn('TXT_realign_intervals', target_creator.tool_outputs)

        indel_realigner = GATKIndelRealigner()
        indel_realigner.add_input_files({
            'BAM': [TestGATK.FILE_BAM],
            'FASTA_REF': [TestGATK.FILE_FASTA_REF],
            'TXT_realign_intervals': target_creator.tool_outputs['TXT_realign_intervals']
        })
        indel_realigner.run(self.__create_dir('indel_realigner'))
        self.verify_output_files(indel_realigner, 'BAM')

    def test_base_quality_score_recalibration(self) -> None:
        """
        Tests the base quality score recalibration: base recalibrator, print reads and analyze covariates.
        :return: None
        """
        table_before = self.__run_base_recalibrator(TestGATK.FILE_BAM, 'recalibrator_before')
        print_reads = GATKPrintReads()
        print_reads.add_input_files({
            'BAM': [TestGATK.FILE_BAM],
            'FASTA_REF': [TestGATK.FILE_FASTA_REF],
            'BQSR': [table_before]
        })
        print_reads.run(self.__create_dir('print_reads'))
        self.verify_output_files(print_reads, 'BAM')

        table_after = self.__run_base_recalibrator(print_reads.tool_outputs['BAM'][0], 'recalibrator_after')
        analyze_covariates = GATKAnalyzeCovariates()
        analyze_covariates.add_input_files({
            'TXT_TABLE_BEFORE': [table_before],
            'TXT_TABLE_AFTER': [table_after],
            'FASTA_REF': [TestGATK.FILE_FASTA_REF]
        })
        # The PDF output is disabled, as it requires R packages (gplots, gsalib) that are not always available
        analyze_covariates.update_parameters(pdf_output=False)
        analyze_covariates.run(self.__create_dir('analyze_covariates'))
        self.assertNotIn('PDF', analyze_covariates.tool_outputs)
        self.verify_output_files(analyze_covariates, 'CSV')

    def test_base_recalibrator_no_known_sites(self) -> None:
        """
        Tests that the base recalibrator fails without known sites.
        :return: None
        """
        base_recalibrator = GATKBaseRecalibrator()
        base_recalibrator.add_input_files({'BAM': [TestGATK.FILE_BAM], 'FASTA_REF': [TestGATK.FILE_FASTA_REF]})
        with self.assertRaises(InvalidToolInputError):
            base_recalibrator.run(self.running_dir)

    def test_haplotype_caller(self) -> None:
        """
        Tests the haplotype caller.
        :return: None
        """
        haplotype_caller = GATKHaplotypeCaller()
        haplotype_caller.add_input_files({'BAM': [TestGATK.FILE_BAM], 'FASTA_REF': [TestGATK.FILE_FASTA_REF]})
        haplotype_caller.update_parameters(bam_output='haplotypes.bam')
        haplotype_caller.run(self.running_dir)
        self.verify_output_files(haplotype_caller, 'VCF')
        self.verify_output_files(haplotype_caller, 'BAM')

    def test_haplotype_caller_genotype_gvcfs(self) -> None:
        """
        Tests the haplotype caller in gVCF mode, followed by genotyping the gVCF.
        :return: None
        """
        haplotype_caller = GATKHaplotypeCaller()
        haplotype_caller.add_input_files({'BAM': [TestGATK.FILE_BAM], 'FASTA_REF': [TestGATK.FILE_FASTA_REF]})
        haplotype_caller.update_parameters(emit_ref_confidence='GVCF', output='sample.g.vcf')
        haplotype_caller.run(self.__create_dir('haplotype_caller'))

        genotype_gvcfs = GATKGenotypeGVCFs()
        genotype_gvcfs.add_input_files({
            'gVCF': haplotype_caller.tool_outputs['VCF'],
            'FASTA_REF': [TestGATK.FILE_FASTA_REF]
        })
        genotype_gvcfs.run(self.__create_dir('genotype_gvcfs'))
        self.verify_output_files(genotype_gvcfs, 'VCF_MultipleSample')

    def test_mutect2(self) -> None:
        """
        Tests MuTect2 in tumor-only mode.
        :return: None
        """
        mutect2 = GATKMuTect2()
        mutect2.add_input_files({'BAM_TUMOR': [TestGATK.FILE_BAM], 'FASTA_REF': [TestGATK.FILE_FASTA_REF]})
        mutect2.update_parameters(output_bam=True)
        mutect2.run(self.running_dir)
        self.verify_output_files(mutect2, 'VCF')
        self.verify_output_files(mutect2, 'BAM')
        self.verify_output_files(mutect2, 'BAI')

    def test_select_variants(self) -> None:
        """
        Tests select variants.
        :return: None
        """
        select_variants = GATKSelectVariants()
        select_variants.add_input_files({'VCF': [TestGATK.FILE_VCF], 'FASTA_REF': [TestGATK.FILE_FASTA_REF]})
        select_variants.update_parameters(select='QUAL>30', select_type_to_include='SNP,MNP')
        select_variants.run(self.running_dir)
        self.verify_output_files(select_variants, 'VCF')

    def test_select_variants_whitespace(self) -> None:
        """
        Tests that select variants fails on a selection expression with whitespace.
        :return: None
        """
        select_variants = GATKSelectVariants()
        select_variants.add_input_files({'VCF': [TestGATK.FILE_VCF], 'FASTA_REF': [TestGATK.FILE_FASTA_REF]})
        select_variants.update_parameters(select='QUAL > 30')
        with self.assertRaises(InvalidParameterError):
            select_variants.run(self.running_dir)

    def test_fasta_alternate_reference_maker(self) -> None:
        """
        Tests the FASTA alternate reference maker, concatenating the extracted intervals.
        :return: None
        """
        # Unsorted on purpose: GATK sorts the intervals
        path_intervals = self.running_dir / 'regions.intervals'
        path_intervals.write_text('phiX174:1401-1450\nphiX174:1-100\n')

        reference_maker = GATKFastaAlternateReferenceMaker()
        reference_maker.add_input_files({
            'VCF': [TestGATK.FILE_VCF],
            'FASTA_REF': [TestGATK.FILE_FASTA_REF],
            'TXT_intervals': [ToolIOFile(path_intervals)]
        })
        reference_maker.update_parameters(concatenate_sequence_segments=True)
        reference_maker.run(self.running_dir)
        self.verify_output_files(reference_maker, 'FASTA')
        self.verify_output_files(reference_maker, 'FASTA_concatenated')

        records_extracted = fastautils.read_as_dict(reference_maker.tool_outputs['FASTA'][0].path)
        self.assertEqual(['phiX174:1-100', 'phiX174:1401-1450'], list(records_extracted.keys()))
        records_concatenated = fastautils.read_as_dict(reference_maker.tool_outputs['FASTA_concatenated'][0].path)
        self.assertEqual(['phiX174'], list(records_concatenated.keys()))
        seq_concatenated = str(records_concatenated['phiX174'].seq)
        self.assertEqual(1450, len(seq_concatenated))
        self.assertEqual('N' * 1300, seq_concatenated[100:1400])
        # The SNP (phiX174:1424 C>A) is applied
        self.assertEqual('A', seq_concatenated[1423])


if __name__ == '__main__':
    unittest.main()
