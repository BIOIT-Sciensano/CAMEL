import unittest

from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.scriptutils.basepipe import basepipeutils
from camel.scripts.enterococcuspipeline import CONFIG_DATA
from camel.scripts.enterococcuspipeline.mainenterococcuspipeline import main
from camel.tests import longRunningTest

CUSTOM_ANALYSES = basepipeutils.get_custom_analyses(CONFIG_DATA)


class TestEnterococcusPipeline(CamelTestSuite):
    """
    Tests for the Enterococcus pipeline.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir()
    input_faecalis_fastq_pe = [
        test_file_dir / 'pipelines' / 'Enterococcus_faecalis-SRR12362697-ds_1.fastq.gz',
        test_file_dir / 'pipelines' / 'Enterococcus_faecalis-SRR12362697-ds_2.fastq.gz'
    ]
    input_faecium_fastq_pe = [
        test_file_dir / 'pipelines' / 'Enterococcus_faecium-SRR12388968-ds_1.fastq.gz',
        test_file_dir / 'pipelines' / 'Enterococcus_faecium-SRR12388968-ds_2.fastq.gz'
    ]
    input_gallinarum_fastq_pe = [
        test_file_dir / 'pipelines' / 'Enterococcus_gallinarum-SRR16344675-ds_1.fastq.gz',
        test_file_dir / 'pipelines' / 'Enterococcus_gallinarum-SRR16344675-ds_2.fastq.gz'
    ]
    input_faecalis_fastq_se = test_file_dir / 'pipelines' / 'Enterococcus_faecalis-SRR17731943_ont-ds.fastq.gz'
    input_faecium_fastq_se = test_file_dir / 'pipelines' / 'Enterococcus_faecium-SRR24726895_ont-ds.fastq.gz'
    input_gallinarum_fastq_se = test_file_dir / 'pipelines' / 'Enterococcus_gallinarum-SRR17662735_ont-ds.fastq.gz'
    input_faecalis_fasta = test_file_dir / 'pipelines' / 'Enterococcus_faecalis-SRR12362697-ds.fasta'
    input_faecium_fasta = test_file_dir / 'pipelines' / 'Enterococcus_faecium-SRR12388968-ds.fasta'
    input_gallinarum_fasta = test_file_dir / 'pipelines' / 'Enterococcus_gallinarum-SRR16344675-ds.fasta'

    @longRunningTest()
    def test_faecalis_blast(self) -> None:
        """
        Tests the Enterococcus pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestEnterococcusPipeline.input_faecalis_fastq_pe[0]),
            str(TestEnterococcusPipeline.input_faecalis_fastq_pe[1]),
            '--input-type', 'illumina',
            '--species', 'faecalis',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_faecalis_kma(self) -> None:
        """
        Tests the Enterococcus pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestEnterococcusPipeline.input_faecalis_fastq_pe[0]),
            str(TestEnterococcusPipeline.input_faecalis_fastq_pe[1]),
            '--input-type', 'illumina',
            '--species', 'faecalis',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_faecium_blast(self) -> None:
        """
        Tests the Enterococcus pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestEnterococcusPipeline.input_faecium_fastq_pe[0]),
            str(TestEnterococcusPipeline.input_faecium_fastq_pe[1]),
            '--input-type', 'illumina',
            '--species', 'faecium',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_faecium_kma(self) -> None:
        """
        Tests the Enterococcus pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestEnterococcusPipeline.input_faecium_fastq_pe[0]),
            str(TestEnterococcusPipeline.input_faecium_fastq_pe[1]),
            '--input-type', 'illumina',
            '--species', 'faecium',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_spp_blast(self) -> None:
        """
        Tests the Enterococcus pipeline for generic Enterococcus with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestEnterococcusPipeline.input_gallinarum_fastq_pe[0]),
            str(TestEnterococcusPipeline.input_gallinarum_fastq_pe[1]),
            '--input-type', 'illumina',
            '--species', 'spp',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_faecalis_fasta(self) -> None:
        """
        Tests the Enterococcus pipeline using FASTA as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fasta', str(TestEnterococcusPipeline.input_faecalis_fasta),
            '--input-type', 'fasta',
            '--species', 'faecalis',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_faecium_fasta(self) -> None:
        """
        Tests the Enterococcus pipeline using FASTA as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fasta', str(TestEnterococcusPipeline.input_faecium_fasta),
            '--input-type', 'fasta',
            '--species', 'faecium',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_spp_fasta(self) -> None:
        """
        Tests the Enterococcus pipeline for generic Enterococcus using FASTA as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fasta', str(TestEnterococcusPipeline.input_gallinarum_fasta),
            '--input-type', 'fasta',
            '--species', 'spp',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_faecalis_ont(self) -> None:
        """
        Tests the Enterococcus pipeline using ONT as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestEnterococcusPipeline.input_faecalis_fastq_se),
            '--input-type', 'ont',
            '--species', 'faecalis',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_faecium_ont(self) -> None:
        """
        Tests the Enterococcus pipeline using ONT as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestEnterococcusPipeline.input_faecium_fastq_se),
            '--input-type', 'ont',
            '--species', 'faecium',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_spp_ont(self) -> None:
        """
        Tests the Enterococcus pipeline for generic Enterococcus using ONT as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestEnterococcusPipeline.input_gallinarum_fastq_se),
            '--input-type', 'ont',
            '--species', 'spp',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_faecalis_kma_ont(self) -> None:
        """
        Tests the Enterococcus pipeline using ONT as input with all assays except for cgMLST and kma detection method
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestEnterococcusPipeline.input_faecalis_fastq_se),
            '--input-type', 'ont',
            '--species', 'faecalis',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_faecium_kma_ont(self) -> None:
        """
        Tests the Enterococcus pipeline using ONT as input with all assays except for cgMLST and kma detection method
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestEnterococcusPipeline.input_faecium_fastq_se),
            '--input-type', 'ont',
            '--species', 'faecium',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_spp_kma_ont(self) -> None:
        """
        Tests the Enterococcus pipeline using ONT as input with all assays except for cgMLST and kma detection method
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestEnterococcusPipeline.input_gallinarum_fastq_se),
            '--input-type', 'ont',
            '--species', 'spp',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
