import unittest

from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.scriptutils.basepipe import basepipeutils
from camel.scripts.salmonellapipeline import CONFIG_DATA
from camel.scripts.salmonellapipeline.mainsalmonellapipeline import (
    main,
)
from camel.tests import longRunningTest

CUSTOM_ANALYSES = basepipeutils.get_custom_analyses(CONFIG_DATA)


class TestSalmonellaPipeline(CamelTestSuite):
    """
    Tests for the Salmonella pipeline.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir()
    input_fastq_illumina_pe = [
        test_file_dir / 'pipelines' / "Salmonella-MB6391-ds_1.fastq.gz",
        test_file_dir / 'pipelines' / "Salmonella-MB6391-ds_2.fastq.gz"
    ]
    input_fastq_se = test_file_dir / 'pipelines' / 'Salmonella-S23BD05337-RBK_ont-ds.fastq.gz'
    input_fasta = test_file_dir / 'pipelines' / "Salmonella-MB6391-ds.fasta"

    @longRunningTest()
    def test_blast_illumina(self) -> None:
        """
        Tests the Salmonella pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        path_fasta_out = self.running_dir / 'out' / 'assembly.fasta'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestSalmonellaPipeline.input_fastq_illumina_pe[0]),
            str(TestSalmonellaPipeline.input_fastq_illumina_pe[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-fasta', str(path_fasta_out),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if 'cgmlst' not in a),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)
        self.assertGreater(path_fasta_out.stat().st_size, 0)

    @longRunningTest()
    def test_blast_illumina_with_json(self) -> None:
        """
        Tests the Salmonella pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        path_json_out = self.running_dir / 'out' / 'output.json'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestSalmonellaPipeline.input_fastq_illumina_pe[0]),
            str(TestSalmonellaPipeline.input_fastq_illumina_pe[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-json', str(path_json_out),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if 'cgmlst' not in a),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)
        self.assertGreater(path_json_out.stat().st_size, 0)

    @longRunningTest()
    def test_kma_illumina(self) -> None:
        """
        Tests the Salmonella pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestSalmonellaPipeline.input_fastq_illumina_pe[0]),
            str(TestSalmonellaPipeline.input_fastq_illumina_pe[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if 'cgmlst' not in a),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_fasta(self) -> None:
        """
        Tests the Salmonella pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fasta', str(TestSalmonellaPipeline.input_fasta),
            '--input-type', 'fasta',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if 'cgmlst' not in a),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_ont(self) -> None:
        """
        Tests the Salmonella pipeline with all assays except for cgMLST , ONT input
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestSalmonellaPipeline.input_fastq_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if 'cgmlst' not in a),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_kma_ont(self) -> None:
        """
        Tests the Salmonella pipeline with all assays except for cgMLST, using KMA and ONT input
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestSalmonellaPipeline.input_fastq_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if 'cgmlst' not in a),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
