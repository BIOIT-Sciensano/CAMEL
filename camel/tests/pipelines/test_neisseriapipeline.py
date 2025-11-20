import unittest
from pathlib import Path

from camel.app.cli import cliutils
from camel.app.core import cameltesthelper
from camel.app.tools.pipelines.sequence_typing.typingdbloader import TypingDBLoader

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.scripts.neisseriapipeline import CONFIG_DATA
from camel.scripts.neisseriapipeline.mainneisseriapipeline import main, CUSTOM_ANALYSES
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
        data_typing = cameltesthelper.extract_from_yaml(CONFIG_DATA, 'sequence_typing')
        for key, scheme_data in data_typing['dbs'].items():
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
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestNeisseriaPipeline.input_fastq_pe[0]),
            str(TestNeisseriaPipeline.input_fastq_pe[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if a != 'cgmlst'),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_neisseria_pipeline_blast_fastp(self) -> None:
        """
        Tests the Neisseria pipeline with fastp trimming and all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestNeisseriaPipeline.input_fastq_pe[0]),
            str(TestNeisseriaPipeline.input_fastq_pe[1]),
            '--input-type', 'illumina',
            '--trimming-method', 'fastp',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if a != 'cgmlst'),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_neisseria_pipeline_blast_with_downsampling(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestNeisseriaPipeline.input_fastq_pe[0]),
            str(TestNeisseriaPipeline.input_fastq_pe[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--cov-max', '5',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if a != 'cgmlst'),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_neisseria_pipeline_kma(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestNeisseriaPipeline.input_fastq_pe[0]),
            str(TestNeisseriaPipeline.input_fastq_pe[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if a != 'cgmlst'),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_neisseria_pipeline_fasta(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fasta', str(TestNeisseriaPipeline.input_fasta),
            '--input-type', 'fasta',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if a != 'cgmlst'),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_neisseria_pipeline_ont(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST with ONT input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestNeisseriaPipeline.input_fastq_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if a != 'cgmlst'),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_neisseria_pipeline_kma_ont(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST, KMA with ONT input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestNeisseriaPipeline.input_fastq_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if a != 'cgmlst'),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
