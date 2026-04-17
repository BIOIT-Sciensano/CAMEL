import unittest

from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.scriptutils.basepipe import basepipeutils
from camel.scripts.yersiniapipeline import CONFIG_DATA
from camel.scripts.yersiniapipeline.mainyersiniapipeline import main
from camel.tests import longRunningTest

CUSTOM_ANALYSES = basepipeutils.get_custom_analyses(CONFIG_DATA)


class TestYersiniaPipeline(CamelTestSuite):
    """
    Tests for the Yersinia pipeline.
    """

    running_dir = None

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines')
    input_enterocolitica_fastq_pe = [
        test_file_dir / 'Yersinia-enterocolitica-S23BD07911_NG_A0183-ds_1.fastq.gz',
        test_file_dir / 'Yersinia-enterocolitica-S23BD07911_NG_A0183-ds_2.fastq.gz'
    ]
    input_pseudotuberculosis_fastq_pe = [
        test_file_dir / 'Yersinia_pseudotuberculosis-S23BD09896_NG_A0586-ds_1.fastq.gz',
        test_file_dir / 'Yersinia_pseudotuberculosis-S23BD09896_NG_A0586-ds_2.fastq.gz'
    ]
    input_enterocolitica_fastq_se = test_file_dir / 'Yersinia-enterocolitica-FAZ88297_ont-ds.fastq.gz'
    input_pseudotuberculosis_fastq_se = test_file_dir / 'Yersinia_pseudotuberculosis-FAZ88297_ont-ds.fastq.gz'
    input_enterocolitica_fasta = test_file_dir / 'Yersinia-enterocolitica-S23BD07911_NG_A0183-ds.fasta'
    input_pseudotuberculosis_fasta = test_file_dir / 'Yersinia_pseudotuberculosis-S23BD09896_NG_A0586-ds.fasta'

    @longRunningTest()
    def test_illumina_blast_with_downsampling(self) -> None:
        """
        Tests the Yersinia pipeline with all assays, except for cgMLST, with downsampling.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestYersiniaPipeline.input_enterocolitica_fastq_pe[0]),
            str(TestYersiniaPipeline.input_enterocolitica_fastq_pe[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--cov-max', '5',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--threads', '4',
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_illumina_enterocolitica_blast(self) -> None:
        """
        Tests the Yersinia pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestYersiniaPipeline.input_enterocolitica_fastq_pe[0]),
            str(TestYersiniaPipeline.input_enterocolitica_fastq_pe[1]),
            '--input-type', 'illumina',
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
    def test_illumina_enterocolitica_kma(self) -> None:
        """
        Tests the Neisseria pipeline with all assays, except for cgMLST,
        with the kma detection method and the TruSeq2 library.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestYersiniaPipeline.input_enterocolitica_fastq_pe[0]),
            str(TestYersiniaPipeline.input_enterocolitica_fastq_pe[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst'))
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_illumina_pseudotuberculosis_blast(self) -> None:
        """
        Tests the Yersinia pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestYersiniaPipeline.input_pseudotuberculosis_fastq_pe[0]),
            str(TestYersiniaPipeline.input_pseudotuberculosis_fastq_pe[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst'))
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_illumina_pseudotuberculosis_kma(self) -> None:
        """
        Tests the Yersinia pipeline with all assays, except for cgMLST,
        with the kma detection method and the TruSeq2 library.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestYersiniaPipeline.input_pseudotuberculosis_fastq_pe[0]),
            str(TestYersiniaPipeline.input_pseudotuberculosis_fastq_pe[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst'))
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_fasta_enterocolitica(self) -> None:
        """
        Tests the Yersinia pipeline with all assays except for cgMLST using FASTA as input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fasta', str(TestYersiniaPipeline.input_enterocolitica_fasta),
            '--input-type', 'fasta',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--kraken2-small-db',
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_fasta_pseudotuberculosis(self) -> None:
        """
        Tests the Yersinia pipeline with all assays except for cgMLST using FASTA as input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fasta', str(TestYersiniaPipeline.input_pseudotuberculosis_fasta),
            '--input-type', 'fasta',
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
    def test_ont_enterocolitica(self) -> None:
        """
        Tests the Yersinia pipeline with all assays except for cgMLST using ONT as input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestYersiniaPipeline.input_enterocolitica_fastq_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst'))
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_ont_pseudotuberculosis(self) -> None:
        """
        Tests the Yersinia pipeline with all assays except for cgMLST using ONT as input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestYersiniaPipeline.input_pseudotuberculosis_fastq_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst'))
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_ont_enterocolitica_kma(self) -> None:
        """
        Tests the Yersinia pipeline with all assays except for cgMLST using ONT as input and kma as detection
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestYersiniaPipeline.input_enterocolitica_fastq_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst'))
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_ont_pseudotuberculosis_kma(self) -> None:
        """
        Tests the Yersinia pipeline with all assays except for cgMLST using ONT as input and kma dectetion
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestYersiniaPipeline.input_pseudotuberculosis_fastq_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst'))
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
