import unittest
from pathlib import Path

from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.app.tools.pipelines.sequence_typing.typingdbloader import TypingDBLoader
from camel.scripts.bacilluspipeline.mainbacilluspipeline import CUSTOM_ANALYSES, main
from camel.tests import longRunningTest


class TestBacillusPipeline(CamelTestSuite):
    """
    Tests for the Bacillus pipeline.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines')
    input_fastq_pe_cereus = [
        test_file_dir / 'Bcereus-SRR7067969-ds_1.fastq.gz',
        test_file_dir / 'Bcereus-SRR7067969-ds_2.fastq.gz'
    ]
    input_fastq_pe_subtilis = [
        test_file_dir / 'Bsubtilis-SRR10568181_1.fastq.gz',
        test_file_dir / 'Bsubtilis-SRR10568181_2.fastq.gz'
    ]
    input_fastq_se_cereus = test_file_dir / 'Bcereus-DRR206405_ont-ds.fastq.gz'
    input_fastq_se_subtilis = test_file_dir / 'Bsubtilis-SRR23725160_ont-ds.fastq.gz'
    input_fasta_subtilis = test_file_dir / 'Bsubtilis-SRR10260289.fasta'
    input_fasta_cereus = test_file_dir / 'Bcereus-D12.fasta'

    def test_bacillus_pipeline_typing_db(self) -> None:
        """
        Checks if the databases for the sequence typing are available.
        :return: None
        """
        sequence_typing_dict = {
             'mlst_cereus': {'path': '/db/sequence_typing/bacillus_cereus/mlst'},
             'mlst_subtilis': {'path': '/db/sequence_typing/bacillus_subtilis/mlst'}}
        for key, scheme_data in sequence_typing_dict.items():
            # Check if scheme exists
            self.assertGreater(Path(scheme_data['path']).stat().st_size, 0)

            # Check if metadata can be loaded
            manager = TypingDBLoader()
            manager.add_input_files({'DIR': [ToolIODirectory(Path(scheme_data['path']))]})
            manager.run(self.running_dir)
            self.assertGreater(len(manager.informs), 0)

    @longRunningTest()
    def test_bacillus_pipeline_subtilis_blast_illumina(self) -> None:
        """
        Tests the Bacillus pipeline with blast-based detection.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'
        tested_analyses = CUSTOM_ANALYSES['common'] + CUSTOM_ANALYSES['subtilis']
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestBacillusPipeline.input_fastq_pe_subtilis[0]),
            str(TestBacillusPipeline.input_fastq_pe_subtilis[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8',
            '--species', 'subtilis'
            '--analyses', ','.join(a for a in tested_analyses if 'cgmlst' not in a),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_bacillus_pipeline_cereus_blast_illumina(self) -> None:
        """
        Tests the Bacillus pipeline with blast-based detection.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'
        tested_analyses = CUSTOM_ANALYSES['common'] + CUSTOM_ANALYSES['cereus']
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestBacillusPipeline.input_fastq_pe_cereus[0]),
            str(TestBacillusPipeline.input_fastq_pe_cereus[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8',
            '--species', 'cereus',
            '--mobsuite-contig-report',
            '--analyses', ','.join(a for a in tested_analyses if 'cgmlst' not in a),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_bacillus_pipeline_subtilis_blast_ont(self) -> None:
        """
        Tests the Bacillus pipeline with blast based detection and ONT data.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'
        tested_analyses = CUSTOM_ANALYSES['common'] + CUSTOM_ANALYSES['subtilis']
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestBacillusPipeline.input_fastq_se_subtilis),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8',
            '--species', 'subtilis',
            '--analyses', ','.join(a for a in tested_analyses if 'cgmlst' not in a),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_bacillus_pipeline_cereus_blast_ont(self) -> None:
        """
        Tests the Bacillus pipeline with blast-based detection and ONT data.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'
        tested_analyses = CUSTOM_ANALYSES['common'] + CUSTOM_ANALYSES['cereus']
        result = cliutils.invoke(main, [
             '--fastq-se', str(TestBacillusPipeline.input_fastq_se_cereus),
             '--input-type', 'ont',
             '--output-html', str(path_report_out),
             '--output-dir', str(path_report_out.parent),
             '--output-tsv', str(path_summary_out),
             '--working-dir', str(self.running_dir),
             '--detection-method', 'blast',
             '--threads', '8',
             '--species', 'cereus',
             '--analyses', ','.join(a for a in tested_analyses if 'cgmlst' not in a),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_bacillus_pipeline_subtilis_blast_hybrid(self) -> None:
        """
        Tests the Bacillus pipeline with blast-based detection and ONT data.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'
        tested_analyses = CUSTOM_ANALYSES['common'] + CUSTOM_ANALYSES['subtilis']
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestBacillusPipeline.test_file_dir / 'Bsubtilis-SRR10260288_50x.fastq.gz'),
            '--fastq-pe',
            str(TestBacillusPipeline.test_file_dir / 'Bsubtilis-SRR10260289_5x_1.fastq.gz'),
            str(TestBacillusPipeline.test_file_dir / 'Bsubtilis-SRR10260289_5x_2.fastq.gz'),
            '--input-type', 'hybrid',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8',
            '--species', 'subtilis',
            '--analyses', ','.join(a for a in tested_analyses if 'cgmlst' not in a),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_bacillus_pipeline_subtilis_fasta(self) -> None:
        """
        Tests the Bacillus pipeline on B. subtilis with blast-based detection and FASTA data.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'
        tested_analyses = CUSTOM_ANALYSES['common']
        result = cliutils.invoke(main, [
            '--fasta', str(TestBacillusPipeline.test_file_dir / 'Bsubtilis-SRR10260289.fasta'),
            '--input-type', 'fasta',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8',
            '--species', 'subtilis',
            '--analyses', ','.join(a for a in tested_analyses if 'cgmlst' not in a),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_bacillus_pipeline_cereus_fasta(self) -> None:
        """
        Tests the Bacillus pipeline on B. cereus with blast-based detection and FASTA data.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fasta', str(TestBacillusPipeline.input_fasta_cereus),
            '--input-type', 'fasta',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8',
            '--species', 'cereus',
            '--analyses', 'btyper'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
