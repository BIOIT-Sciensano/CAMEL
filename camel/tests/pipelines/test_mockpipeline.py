import unittest
from pathlib import Path

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.scripts.mockpipeline.mainmockpipeline import MainMockPipeline
from camel.tests import longRunningTest


class TestMockPipeline(CamelTestSuite):
    """
    Tests for the mock pipeline.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines', 'mock_pipeline')
    input_ilmn_pe = [test_file_dir / 'ecoli_10k_ilmn_1.fastq.gz', test_file_dir / 'ecoli_10k_ilmn_2.fastq.gz']
    input_ont_se = test_file_dir / 'ecoli_10k_ont.fastq.gz'
    input_fasta = test_file_dir / 'ecoli_10k.fasta'
    input_vcf = test_file_dir / 'variants-ecoli_10k_ilmn-all.vcf'

    def test_mock_pipeline_illumina_trimmomatic(self) -> None:
        """
        Tests the mock pipeline with Illumina input data and trimmomatic.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'
        pipeline = MainMockPipeline([
            '--input-type', 'illumina',
            '--fastq-pe', *[str(x) for x in TestMockPipeline.input_ilmn_pe],
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--trimming-method', 'trimmomatic',
            '--detection-method', 'blast',
            '--ncbi-amr',
            '--human-read-scrubbing',
            '--threads', '8'
        ])
        pipeline.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_mock_pipeline_illumina_fastp(self) -> None:
        """
        Tests the mock pipeline with Illumina input data and fastp trimming.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'
        path_fasta_out = Path(self.running_dir) / 'out' / 'assembly.fasta'
        pipeline = MainMockPipeline([
            '--input-type', 'illumina',
            '--fastq-pe', *[str(x) for x in TestMockPipeline.input_ilmn_pe],
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--output-fasta', str(path_fasta_out),
            '--working-dir', str(self.running_dir),
            '--trimming-method', 'fastp',
            '--detection-method', 'blast',
            '--ncbi-amr',
            '--human-read-scrubbing',
            '--threads', '8'
        ])
        pipeline.run()
        self.assertGreater(path_report_out.stat().st_size, 0)
        self.assertGreater(path_fasta_out.stat().st_size, 0)

    def test_mock_pipeline_illumina_fastp_with_json(self) -> None:
        """
        Tests the mock pipeline with Illumina input data and fastp trimming and JSON export.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'
        path_json_out = Path(self.running_dir) / 'out' / 'output.json'
        pipeline = MainMockPipeline([
            '--input-type', 'illumina',
            '--fastq-pe', *[str(x) for x in TestMockPipeline.input_ilmn_pe],
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--output-json', str(path_json_out),
            '--working-dir', str(self.running_dir),
            '--trimming-method', 'fastp',
            '--detection-method', 'blast',
            '--ncbi-amr',
            '--human-read-scrubbing',
            '--threads', '8'
        ])
        pipeline.run()
        self.assertGreater(path_report_out.stat().st_size, 0)
        self.assertGreater(path_json_out.stat().st_size, 0)

    def test_mock_pipeline_illumina_kma(self) -> None:
        """
        Tests the mock pipeline with Illumina input data and KMA detection.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'
        pipeline = MainMockPipeline([
            '--input-type', 'illumina',
            '--fastq-pe', *[str(x) for x in TestMockPipeline.input_ilmn_pe],
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--ncbi-amr',
            '--threads', '8'
        ])
        pipeline.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_mock_pipeline_ont(self) -> None:
        """
        Tests the mock pipeline with ONT input data.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'

        pipeline = MainMockPipeline([
            '--fastq-se', str(TestMockPipeline.input_ont_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--ncbi-amr',
            '--threads', '8',
            '--confindr'
        ])
        pipeline.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_mock_pipeline_ont_with_downsampling(self) -> None:
        """
        Tests the mock pipeline with ONT input data.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'

        pipeline = MainMockPipeline([
            '--fastq-se', str(TestMockPipeline.input_ont_se),
            '--input-type', 'ont',
            '--cov-max', '50',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--ncbi-amr',
            '--threads', '8',
        ])
        pipeline.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_mock_pipeline_ont_with_scrubbing_and_filters_params(self) -> None:
        """
        Tests the mock pipeline with ONT input data and read scrubbing enabled and alternate read filtering parameters.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'

        pipeline = MainMockPipeline([
            '--fastq-se', str(TestMockPipeline.input_ont_se),
            '--input-type', 'ont',
            '--cov-max', '50',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--ncbi-amr',
            '--human-read-scrubbing',
            '--ont-min-qual', '11',
            '--ont-min-length', '750',
            '--threads', '8',
        ])
        pipeline.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_mock_pipeline_hybrid(self) -> None:
        """
        Tests the mock pipeline with hybrid input data (Illumina + ONT).
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'

        pipeline = MainMockPipeline([
            '--fastq-se', str(TestMockPipeline.input_ont_se),
            '--fastq-pe', *[str(x) for x in TestMockPipeline.input_ilmn_pe],
            '--input-type', 'hybrid',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--ncbi-amr',
            '--human-read-scrubbing',
            '--threads', '8',
        ])
        pipeline.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_mock_pipeline_fasta(self) -> None:
        """
        Tests the mock pipeline with FASTA input data.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'

        pipeline = MainMockPipeline([
            '--fasta', str(TestMockPipeline.input_fasta),
            '--input-type', 'fasta',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--ncbi-amr',
            '--snpit',
            '--human-read-scrubbing',
            '--threads', '8',
        ])
        pipeline.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_mock_pipeline_fasta_with_kraken2(self) -> None:
        """
        Tests the mock pipeline with FASTA input data with the Kraken2 analysis enabled.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'

        pipeline = MainMockPipeline([
            '--fasta', str(TestMockPipeline.input_fasta),
            '--input-type', 'fasta',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--kraken2',
            '--ncbi-amr',
            '--human-read-scrubbing',
            '--threads', '8',
        ])
        pipeline.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_mock_pipeline_fasta_with_vcf(self) -> None:
        """
        Tests the mock pipeline with FASTA and VCF input data.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'

        pipeline = MainMockPipeline([
            '--fasta', str(TestMockPipeline.input_fasta),
            '--vcf-unfiltered', str(TestMockPipeline.input_vcf),
            '--input-type', 'fasta_with_vcf',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--ncbi-amr',
            '--snpit',
            '--human-read-scrubbing',
            '--threads', '8',
        ])
        pipeline.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
