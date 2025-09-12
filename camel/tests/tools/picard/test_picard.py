import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.picard.addorreplacereadgroups import AddOrReplaceReadGroups
from camel.app.tools.picard.calculatereadgroupchecksum import CalculateReadGroupChecksum
from camel.app.tools.picard.collectmultiplemetrics import CollectMultipleMetrics
from camel.app.tools.picard.collectqualityyieldmetrics import CollectQualityYieldMetrics
from camel.app.tools.picard.collectrawwgsmetrics import CollectRawWgsMetrics
from camel.app.tools.picard.collectvariantcallingmetrics import CollectVariantCallingMetrics
from camel.app.tools.picard.collectwgsmetrics import CollectWgsMetrics
from camel.app.tools.picard.createsequencedictionary import CreateSequenceDictionary
from camel.app.tools.picard.fastqtosam import FastqToSam
from camel.app.tools.picard.gatherbamfiles import GatherBamFiles
from camel.app.tools.picard.intervallisttools import IntervalListTools
from camel.app.tools.picard.markduplicates import MarkDuplicates
from camel.app.tools.picard.mergebamalignment import MergeBamAlignment
from camel.app.tools.picard.mergevcfs import MergeVCFs
from camel.app.tools.picard.samtofastq import SamToFastq
from camel.app.tools.picard.setnmmdanduqtags import SetNmMdAndUqTags
from camel.app.tools.picard.sortsam import SortSam
from camel.app.tools.picard.validatesamfile import ValidateSamFile
from camel.tests import minOSVersion


class TestPicard(CamelTestSuite):
    """
    Tests the Picard tool suite.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('picard')
    FILE_BAM = ToolIOFile(test_file_dir / 'aln_rg.bam')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'reference.fasta')
    FILE_VCFdb = ToolIOFile(test_file_dir / 'Homo_sapiens_assembly38.dbsnp138_chr22.vcf.gz')
    FILE_BAM_SORTED = ToolIOFile(test_file_dir / 'sorted.bam')
    FILE_VCF = ToolIOFile(test_file_dir / 'unfiltered_variants-myco.vcf')
    FILE_VCF1 = ToolIOFile(test_file_dir / 'var1.vcf')
    FILE_VCF2 = ToolIOFile(test_file_dir / 'var2.vcf')
    FILE_FASTQ_R1 = ToolIOFile(test_file_dir / 'r_1.fq')
    FILE_FASTQ_R2 = ToolIOFile(test_file_dir / 'r_2.fq')
    FILE_BAM1 = ToolIOFile(test_file_dir / 'aln1.bam')
    FILE_BAM2 = ToolIOFile(test_file_dir / 'aln2.bam')
    FILE_uBAM = ToolIOFile(test_file_dir / 'unmapped.bam')
    FILE_VCF_human = ToolIOFile(test_file_dir / 'NA12877_chr22.g.vcf.gz')
    FILE_BAM_human = ToolIOFile(test_file_dir / 'readgroup_updated.bam')
    FILE_REF_human = ToolIOFile(test_file_dir / 'Homo_sapiens_assembly38_chr22.fasta')

    def test_picard_addorreplacereadgroups(self) -> None:
        """
        Test Picard AddOrReplaceReadgroups
        :return: None
        """
        picard_addorreplacereadgroups = AddOrReplaceReadGroups()
        picard_addorreplacereadgroups.add_input_files({
            'BAM': [TestPicard.FILE_BAM]
        })
        picard_addorreplacereadgroups.run(self.running_dir)
        self.verify_output_files(picard_addorreplacereadgroups, 'BAM')

    def test_picard_calculatereadgroupchecksum(self) -> None:
        """
        Test Picard CalculateReadgroupChecksum
        :return: None
        """
        picard_calculatereadgroupchecksum = CalculateReadGroupChecksum()
        picard_calculatereadgroupchecksum.add_input_files({
            'BAM': [TestPicard.FILE_BAM]
        })
        picard_calculatereadgroupchecksum.run(self.running_dir)
        self.verify_output_files(picard_calculatereadgroupchecksum, 'TXT_checksum')

    @minOSVersion('jammy')
    def test_picard_collectmultiplemetrics(self) -> None:
        """
        Test Picard CollectMultipleMetrics
        :return: None
        """
        picard_collectmultiplemetrics = CollectMultipleMetrics()
        picard_collectmultiplemetrics.add_input_files({
            'BAM': [TestPicard.FILE_BAM_human],
            'FASTA_REF': [TestPicard.FILE_REF_human],
        })
        picard_collectmultiplemetrics.update_parameters(
            assume_sorted='true',
            reset_metrics='null',
            metrics_CollectAlignmentSummaryMetrics='CollectAlignmentSummaryMetrics',
            metrics_CollectGcBiasMetrics='CollectGcBiasMetrics',
            metrics_CollectInsertSizeMetrics='CollectInsertSizeMetrics',
            metrics_QualityScoreDistribution='QualityScoreDistribution ',
            metrics_MeanQualityByCycle='MeanQualityByCycle',
            metrics_CollectBaseDistributionByCycle='CollectBaseDistributionByCycle',
            metrics_CollectSequencingArtifactMetrics='CollectSequencingArtifactMetrics',
            metrics_CollectQualityYieldMetrics='CollectQualityYieldMetrics'
        )
        picard_collectmultiplemetrics.run(self.running_dir)

        expected_output = ['TXT_AlignmentSummary', 'TXT_GcBias', 'TXT_GcBiasSummary', 'TXT_GcBiasFigure',
                           'TXT_InsertSize', 'TXT_InsertSizeFigure', 'TXT_QualityDistribution',
                           'TXT_QualityDistributionFigure', 'TXT_QualityByCycle',
                           'TXT_QualityByCycleFigure', 'TXT_BaseDistributionByCycle',
                           'TXT_BaseDistributionByCycleFigure', 'TXT_SequencingArtefactDetail',
                           'TXT_SequencingArtefactSummary', 'TXT_SequencingArtefactErrorSummary',
                           'TXT_SequencingArtefactPreAdapterDetail', 'TXT_SequencingArtefactPreAdapterSummary',
                           'TXT_QualityYield']

        for expected_file in expected_output:
            self.verify_output_files(picard_collectmultiplemetrics, expected_file)

    def test_picard_collectqualityyieldmetrics(self) -> None:
        """
        Test Picard CollectQualityYieldMetrics
        :return: None
        """
        picard_collectqualityyieldmetrics = CollectQualityYieldMetrics()
        picard_collectqualityyieldmetrics.add_input_files({
            'BAM': [TestPicard.FILE_BAM]
        })
        picard_collectqualityyieldmetrics.run(self.running_dir)
        self.verify_output_files(picard_collectqualityyieldmetrics, 'TXT')

    def test_picard_collectrawwgsmetrics(self) -> None:
        """
        Test Picard CollectRawWgsMetrics
        :return: None
        """
        picard_collectrawwgsmetrics = CollectRawWgsMetrics()
        picard_collectrawwgsmetrics.add_input_files({
            'BAM': [TestPicard.FILE_BAM],
            'FASTA_REF': [TestPicard.FILE_FASTA_REF]
        })
        picard_collectrawwgsmetrics.update_parameters(sample_size=10)
        picard_collectrawwgsmetrics.run(self.running_dir)
        self.verify_output_files(picard_collectrawwgsmetrics, 'TXT_metrics')

    def test_picard_collectwgsmetrics(self) -> None:
        """
        Test Picard CollectWgsMetrics
        :return: None
        """
        picard_collectwgsmetrics = CollectWgsMetrics()
        picard_collectwgsmetrics.add_input_files({
            'BAM': [TestPicard.FILE_BAM],
            'FASTA_REF': [TestPicard.FILE_FASTA_REF]
        })
        picard_collectwgsmetrics.update_parameters(sample_size=10)
        picard_collectwgsmetrics.run(self.running_dir)
        self.verify_output_files(picard_collectwgsmetrics, 'TXT_metrics')

    def test_picard_collectvariantcallingmetrics(self) -> None:
        """
        Test Picard CollectVariantCallingMetrics
        :return: None
        """
        picard_collectvariantcallingmetrics = CollectVariantCallingMetrics()
        picard_collectvariantcallingmetrics.add_input_files({
            'VCF': [TestPicard.FILE_VCF_human],
            'VCF_dbsnp': [TestPicard.FILE_VCFdb]
        })
        picard_collectvariantcallingmetrics.run(self.running_dir)
        self.verify_output_files(picard_collectvariantcallingmetrics, 'TXT_report', 2)

    def test_picard_createsequencedictionary(self) -> None:
        """
        Test Picard CreateSequenceDictionary
        :return: None
        """
        picard_createsequencedictionary = CreateSequenceDictionary()
        picard_createsequencedictionary.add_input_files({
            'FASTA_REF': [TestPicard.FILE_FASTA_REF]
        })
        picard_createsequencedictionary.run(self.running_dir)
        self.verify_output_files(picard_createsequencedictionary, 'FASTA_REF')

    def test_picard_fastqtosam(self) -> None:
        """
        Test Picard FastqToSam
        :return: None
        """
        picard_fastqtosam = FastqToSam()
        picard_fastqtosam.add_input_files({
            'FASTQ_PE': [TestPicard.FILE_FASTQ_R1, TestPicard.FILE_FASTQ_R2]
        })
        picard_fastqtosam.run(self.running_dir)
        self.verify_output_files(picard_fastqtosam, 'BAM')

    def test_picard_gatherbamfiles(self) -> None:
        """
        Test Picard GatherBamFiles
        :return: None
        """
        picard_gatherbamfiles = GatherBamFiles()
        picard_gatherbamfiles.add_input_files({
            'BAMs': [TestPicard.FILE_BAM1, TestPicard.FILE_BAM2]
        })
        picard_gatherbamfiles.run(self.running_dir)
        self.verify_output_files(picard_gatherbamfiles, 'BAM')

    def test_picard_intervallisttools(self) -> None:
        """
        Test Picard IntervalListTools
        :return: None
        """
        picard_intervallisttools = IntervalListTools()
        picard_intervallisttools.add_input_files({
            'VCF': [TestPicard.FILE_VCF]
        })

        picard_intervallisttools.run(self.running_dir)
        self.verify_output_files(picard_intervallisttools, 'TXT_intervalLists')

    def test_picard_markduplicates(self) -> None:
        """
        Test Picard MarkDuplicates
        :return: None
        """
        picard_markduplicates = MarkDuplicates()
        picard_markduplicates.add_input_files({
            'BAM': [TestPicard.FILE_BAM_SORTED],
            'FASTA_REF': [TestPicard.FILE_FASTA_REF]
        })
        picard_markduplicates.run(self.running_dir)
        self.verify_output_files(picard_markduplicates, 'BAM')
        self.verify_output_files(picard_markduplicates, 'METRICS')

    def test_picard_mergebamalignment(self) -> None:
        """
        Test Picard MergeBamAlignment
        :return: None
        """
        picard_mergebamalignment = MergeBamAlignment()
        picard_mergebamalignment.add_input_files({
            'BAM_UNMAPPED': [TestPicard.FILE_uBAM],
            'BAM_ALIGNED': [TestPicard.FILE_BAM],
            'FASTA_REF': [TestPicard.FILE_FASTA_REF]
        })
        picard_mergebamalignment.update_parameters(
            attributes_to_remove_multi="NM,MD"
        )
        picard_mergebamalignment.run(self.running_dir)
        self.verify_output_files(picard_mergebamalignment, 'BAM')

    def test_picard_mergevcfs(self) -> None:
        """
        Test Picard MergeVCFs
        :return: None
        """
        picard_mergevcfs = MergeVCFs()
        picard_mergevcfs.add_input_files({
            'VCF': [TestPicard.FILE_VCF1, TestPicard.FILE_VCF2],
        })
        picard_mergevcfs.run(self.running_dir)
        self.verify_output_files(picard_mergevcfs, 'VCF')

    def test_picard_samtofastq(self) -> None:
        """
        Test Picard SamToFastq
        :return: None
        """
        picard_samtofastq = SamToFastq()
        picard_samtofastq.add_input_files({
            'BAM': [TestPicard.FILE_BAM],
        })
        picard_samtofastq.run(self.running_dir)
        self.verify_output_files(picard_samtofastq, 'FASTQ', 2)

    def test_picard_setnmmdanduqtags(self) -> None:
        """
        Test Picard SetNmMdAndUqTags
        :return: None
        """
        picard_setnmmdanduqtags = SetNmMdAndUqTags()
        picard_setnmmdanduqtags.add_input_files({
            'BAM': [TestPicard.FILE_BAM_SORTED],
            'FASTA_REF': [TestPicard.FILE_FASTA_REF]
        })
        picard_setnmmdanduqtags.run(self.running_dir)
        self.verify_output_files(picard_setnmmdanduqtags, 'BAM')

    def test_picard_sortsam(self) -> None:
        """
        Test Picard SortSam
        :return: None
        """
        picard_sortsam = SortSam()
        picard_sortsam.add_input_files({
            'BAM': [TestPicard.FILE_BAM],
            'FASTA_REF': [TestPicard.FILE_FASTA_REF]
        })
        picard_sortsam.run(self.running_dir)
        self.verify_output_files(picard_sortsam, 'BAM')

    def test_picard_validatesamfile(self) -> None:
        """
        Test Picard ValidateSamFile
        :return: None
        """
        picard_validatesamfile = ValidateSamFile()
        picard_validatesamfile.add_input_files({
            'BAM': [TestPicard.FILE_BAM],
            'FASTA_REF': [TestPicard.FILE_FASTA_REF]
        })
        picard_validatesamfile.run(self.running_dir)
        self.verify_output_files(picard_validatesamfile, 'TXT_report')


if __name__ == '__main__':
    unittest.main()
