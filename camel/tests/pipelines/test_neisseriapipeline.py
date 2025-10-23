import unittest
from pathlib import Path
from camel.app.tools.pipelines.sequence_typing.typingdbloader import TypingDBLoader
import yaml

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.scripts.neisseriapipeline import CONFIG_DATA
from camel.scripts.neisseriapipeline.mainneisseriapipeline import MainNeisseriaPipeline
from camel.tests import longRunningTest


class TestNeisseriaPipeline(CamelTestSuite):
    """
    Tests for the Neisseria pipeline.
    """
    running_dir = None

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines')
    input_fastq_pe = [
        test_file_dir / 'Neisseria-2011-006_S6-ds_1.fastq.gz',
        test_file_dir / 'Neisseria-2011-006_S6-ds_2.fastq.gz'
    ]
    input_fastq_se = test_file_dir / 'Neisseria-S16BD06814-RPB_ont-ds.fastq.gz'
    input_fasta = test_file_dir / 'Neisseria-2011-006_S6-ds.fasta'

    def test_neisseria_pipeline_typing_db(self) -> None:
        """
        Checks if the databases for the sequence typing are available.
        :return: None
        """
        with open(CONFIG_DATA) as handle_in:
            config_data = yaml.safe_load(handle_in)

        for key, scheme_data in config_data['sequence_typing']['dbs'].items():
            # Check if scheme exists
            self.assertGreater(Path(scheme_data['path']).stat().st_size, 0)

            # Check if metadata can be loaded
            manager = TypingDBLoader()
            manager.add_input_files({'DIR': [ToolIODirectory(Path(scheme_data['path']))]})
            manager.run(self.running_dir)
            self.assertGreater(len(manager.informs), 0)

    @longRunningTest()
    def test_neisseria_pipeline_blast(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestNeisseriaPipeline.input_fastq_pe[0]), str(TestNeisseriaPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ] + [f"--{a.replace('_', '-')}" for a in MainNeisseriaPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainNeisseriaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_neisseria_pipeline_blast_fastp(self) -> None:
        """
        Tests the Neisseria pipeline with fastp trimming and all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestNeisseriaPipeline.input_fastq_pe[0]), str(TestNeisseriaPipeline.input_fastq_pe[1]),
            '--trimming-method', 'fastp',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ] + [f"--{a.replace('_', '-')}" for a in MainNeisseriaPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainNeisseriaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_neisseria_pipeline_blast_with_downsampling(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestNeisseriaPipeline.input_fastq_pe[0]), str(TestNeisseriaPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--cov-max', '5.0',
        ] + [f"--{a.replace('_', '-')}" for a in MainNeisseriaPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainNeisseriaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_neisseria_pipeline_kma(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestNeisseriaPipeline.input_fastq_pe[0]), str(TestNeisseriaPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--library', 'TruSeq2'
        ] + [f"--{a.replace('_', '-')}" for a in MainNeisseriaPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainNeisseriaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_neisseria_pipeline_fasta(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fasta', str(TestNeisseriaPipeline.input_fasta),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--input-type', 'fasta',
            '--detection-method', 'blast',
        ] + [f"--{a.replace('_', '-')}" for a in MainNeisseriaPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainNeisseriaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_neisseria_pipeline_ont(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST with ONT input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-se', str(TestNeisseriaPipeline.input_fastq_se),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--input-type', 'ont',
            '--detection-method', 'blast',
        ] + [f"--{a.replace('_', '-')}" for a in MainNeisseriaPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainNeisseriaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_neisseria_pipeline_kma_ont(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST, KMA with ONT input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-se', str(TestNeisseriaPipeline.input_fastq_se),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--input-type', 'ont',
            '--detection-method', 'kma',
        ] + [f"--{a.replace('_', '-')}" for a in MainNeisseriaPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainNeisseriaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
