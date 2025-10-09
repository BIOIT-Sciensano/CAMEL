import unittest
from pathlib import Path

import tempfile
import yaml

from camel.app.config import config
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.scripts.yersiniapipeline import CONFIG_DATA
from camel.scripts.yersiniapipeline.mainyersiniapipeline import MainYersiniaPipeline
from camel.tests import longRunningTest


class TestYersiniaPipeline(unittest.TestCase):
    """
    Tests for the Yersinia pipeline.
    """

    running_dir = None

    # Input files
    test_file_dir = Path(config.dir_testdata, 'pipelines')
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

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = Path(tempfile.mkdtemp(None, 'camel_', config.dir_temp))

    def test_yersinia_pipeline_typing_db(self) -> None:
        """
        Checks if the databases for the sequence typing are available.
        :return: None
        """
        from camel.app.tools.pipelines.sequence_typing.typingdbloader import TypingDBLoader
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
    def test_yersinia_pipeline_blast_with_downsampling(self) -> None:
        """
        Tests the Yersinia pipeline with all assays, except for cgMLST, with downsampling.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestYersiniaPipeline.input_enterocolitica_fastq_pe[0]), str(TestYersiniaPipeline.input_enterocolitica_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--cov-max', '5.0',
        ] + [f"--{a.replace('_', '-')}" for a in MainYersiniaPipeline.CUSTOM_ANALYSES if not a.startswith('cgmlst')]
        main = MainYersiniaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_yersinia_pipeline_enterocolitica_blast(self) -> None:
        """
        Tests the Yersinia pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestYersiniaPipeline.input_enterocolitica_fastq_pe[0]), str(TestYersiniaPipeline.input_enterocolitica_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ] + [f"--{a.replace('_', '-')}" for a in MainYersiniaPipeline.CUSTOM_ANALYSES if not a.startswith('cgmlst')]
        main = MainYersiniaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_yersinia_pipeline_enterocolitica_kma(self) -> None:
        """
        Tests the Neisseria pipeline with all assays, except for cgMLST,
        with the kma detection method and the TruSeq2 library.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                '--fastq-pe', str(TestYersiniaPipeline.input_enterocolitica_fastq_pe[0]),
                str(TestYersiniaPipeline.input_enterocolitica_fastq_pe[1]),
                '--output-html', str(path_report_out),
                '--output-dir', str(path_report_out.parent),
                '--output-tsv', str(path_summary_out),
                '--working-dir', str(self.running_dir),
                '--detection-method', 'kma',
                '--library', 'TruSeq2'
            ] + [f"--{a.replace('_', '-')}" for a in MainYersiniaPipeline.CUSTOM_ANALYSES if not a.startswith('cgmlst')]
        main = MainYersiniaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_yersinia_pipeline_pseudotuberculosis_blast(self) -> None:
        """
        Tests the Yersinia pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestYersiniaPipeline.input_pseudotuberculosis_fastq_pe[0]), str(TestYersiniaPipeline.input_pseudotuberculosis_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ] + [f"--{a.replace('_', '-')}" for a in MainYersiniaPipeline.CUSTOM_ANALYSES if not a.startswith('cgmlst')]
        main = MainYersiniaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_yersinia_pipeline_pseudotuberculosis_kma(self) -> None:
        """
        Tests the Yersinia pipeline with all assays, except for cgMLST,
        with the kma detection method and the TruSeq2 library.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                '--fastq-pe', str(TestYersiniaPipeline.input_pseudotuberculosis_fastq_pe[0]),
                str(TestYersiniaPipeline.input_pseudotuberculosis_fastq_pe[1]),
                '--output-html', str(path_report_out),
                '--output-dir', str(path_report_out.parent),
                '--output-tsv', str(path_summary_out),
                '--working-dir', str(self.running_dir),
                '--detection-method', 'kma',
                '--library', 'TruSeq2'
            ] + [f"--{a.replace('_', '-')}" for a in MainYersiniaPipeline.CUSTOM_ANALYSES if not a.startswith('cgmlst')]
        main = MainYersiniaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_yersinia_pipeline_enterocolitica_fasta(self) -> None:
        """
        Tests the Yersinia pipeline with all assays except for cgMLST using FASTA as input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fasta', str(TestYersiniaPipeline.input_enterocolitica_fasta),
                   '--input-type', 'fasta',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir)
               ] + [f"--{a.replace('_', '-')}" for a in MainYersiniaPipeline.CUSTOM_ANALYSES if
                    not a.startswith('cgmlst')]
        main = MainYersiniaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_yersinia_pipeline_pseudotuberculosis_fasta(self) -> None:
        """
        Tests the Yersinia pipeline with all assays except for cgMLST using FASTA as input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fasta', str(TestYersiniaPipeline.input_pseudotuberculosis_fasta),
                   '--input-type', 'fasta',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir)
               ] + [f"--{a.replace('_', '-')}" for a in MainYersiniaPipeline.CUSTOM_ANALYSES if
                    not a.startswith('cgmlst')]
        main = MainYersiniaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_yersinia_pipeline_enterocolitica_ont(self) -> None:
        """
        Tests the Yersinia pipeline with all assays except for cgMLST using ONT as input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fastq-se', str(TestYersiniaPipeline.input_enterocolitica_fastq_se),
                   '--input-type', 'ont',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir)
               ] + [f"--{a.replace('_', '-')}" for a in MainYersiniaPipeline.CUSTOM_ANALYSES if
                    not a.startswith('cgmlst')]
        main = MainYersiniaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_yersinia_pipeline_pseudotuberculosis_ont(self) -> None:
        """
        Tests the Yersinia pipeline with all assays except for cgMLST using ONT as input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fastq-se', str(TestYersiniaPipeline.input_pseudotuberculosis_fastq_se),
                   '--input-type', 'ont',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir)
               ] + [f"--{a.replace('_', '-')}" for a in MainYersiniaPipeline.CUSTOM_ANALYSES if
                    not a.startswith('cgmlst')]
        main = MainYersiniaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_yersinia_pipeline_enterocolitica_kma_ont(self) -> None:
        """
        Tests the Yersinia pipeline with all assays except for cgMLST using ONT as input and kma as detection
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fastq-se', str(TestYersiniaPipeline.input_enterocolitica_fastq_se),
                   '--input-type', 'ont',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir),
                   '--detection-method', 'kma',
               ] + [f"--{a.replace('_', '-')}" for a in MainYersiniaPipeline.CUSTOM_ANALYSES if
                    not a.startswith('cgmlst')]
        main = MainYersiniaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_yersinia_pipeline_pseudotuberculosis_kma_ont(self) -> None:
        """
        Tests the Yersinia pipeline with all assays except for cgMLST using ONT as input and kma dectetion
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fastq-se', str(TestYersiniaPipeline.input_pseudotuberculosis_fastq_se),
                   '--input-type', 'ont',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir),
                   '--detection-method', 'kma',
               ] + [f"--{a.replace('_', '-')}" for a in MainYersiniaPipeline.CUSTOM_ANALYSES if
                    not a.startswith('cgmlst')]
        main = MainYersiniaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
