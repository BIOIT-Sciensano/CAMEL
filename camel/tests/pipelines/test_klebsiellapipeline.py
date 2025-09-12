import unittest
from pathlib import Path

import yaml

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.config import config
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.tools.pipelines.sequence_typing.locussetmanager import LocusSetManager
from camel.scripts.klebsiellapipeline import CONFIG_DATA
from camel.scripts.klebsiellapipeline.mainklebsiellapipeline import MainKlebsiellaPipeline
from camel.tests import longRunningTest


class TestKlebsiellaPipeline(CamelTestSuite):
    """
    Tests for the Klebsiella pipeline.
    """
    running_dir = None

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines')
    input_fastq_pe = [
        test_file_dir / 'Kpneumoniae-SRR4046826-ds_1.fastq.gz',
        test_file_dir / 'Kpneumoniae-SRR4046826-ds_2.fastq.gz'
    ]
    input_fastq_se = test_file_dir / 'Kpneumoniae-FAZ63816_ont-ds.fastq.gz'
    input_fasta = test_file_dir / 'Kpneumoniae-SRR4046826-ds.fasta'

    def test_klebsiella_pipeline_typing_db(self) -> None:
        """
        Checks if the databases for the sequence typing are available.
        :return: None
        """
        with open(CONFIG_DATA) as handle_in:
            config_data = yaml.safe_load(handle_in)

        for key, scheme_data in config_data['sequence_typing'].items():
            # Check if scheme exists
            self.assertGreater(Path(scheme_data['path']).stat().st_size, 0)

            # Check if metadata can be loaded
            manager = LocusSetManager()
            manager.add_input_files({'DIR': [ToolIODirectory(Path(scheme_data['path']))]})
            manager.run(self.running_dir)
            self.assertGreater(len(manager.informs), 0)

    @longRunningTest()
    def test_klebsiella_pipeline_blast(self) -> None:
        """
        Tests the Klebsiella pipeline with blast based detection.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestKlebsiellaPipeline.input_fastq_pe[0]), str(TestKlebsiellaPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8'
        ] + [
            f"--{a.replace('_', '-')}" for a in MainKlebsiellaPipeline.CUSTOM_ANALYSES if a not in (
                'cgmlst', 'scgmlst')]
        main = MainKlebsiellaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_klebsiella_pipeline_kma(self) -> None:
        """
        Tests the Klebsiella pipeline with KMA-based detection.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestKlebsiellaPipeline.input_fastq_pe[0]), str(TestKlebsiellaPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--threads', '8'
        ] + [
            f"--{a.replace('_', '-')}" for a in MainKlebsiellaPipeline.CUSTOM_ANALYSES if a not in (
                'cgmlst', 'scgmlst')]
        main = MainKlebsiellaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_klebsiella_pipeline_fasta(self) -> None:
        """
        Tests the Klebsiella pipeline using FASTA as input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fasta', str(TestKlebsiellaPipeline.input_fasta),
                   '--input-type', 'fasta',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir),
                   '--detection-method', 'blast',
                   '--threads', '8'
               ] + [
                   f"--{a.replace('_', '-')}" for a in MainKlebsiellaPipeline.CUSTOM_ANALYSES if a not in (
                    'cgmlst', 'scgmlst')]
        main = MainKlebsiellaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_klebsiella_pipeline_ont(self) -> None:
        """
        Tests the Klebsiella pipeline using ONT as input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fastq-se', str(TestKlebsiellaPipeline.input_fastq_se),
                   '--input-type', 'ont',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir),
                   '--detection-method', 'blast',
                   '--threads', '8'
               ] + [
                   f"--{a.replace('_', '-')}" for a in MainKlebsiellaPipeline.CUSTOM_ANALYSES if a not in (
                    'cgmlst', 'scgmlst')]
        main = MainKlebsiellaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_klebsiella_pipeline_ont_kma(self) -> None:
        """
        Tests the Klebsiella pipeline using ONT as input and kma as detection method
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
               '--fastq-se', str(TestKlebsiellaPipeline.input_fastq_se),
               '--input-type', 'ont',
               '--output-html', str(path_report_out),
               '--output-dir', str(path_report_out.parent),
               '--output-tsv', str(path_summary_out),
               '--working-dir', str(self.running_dir),
               '--detection-method', 'kma',
               '--threads', '8'
           ] + [
               f"--{a.replace('_', '-')}" for a in MainKlebsiellaPipeline.CUSTOM_ANALYSES if a not in (
                'cgmlst', 'scgmlst')]
        main = MainKlebsiellaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
