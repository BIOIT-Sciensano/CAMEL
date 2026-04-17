import unittest

from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.scriptutils.basepipe import basepipeutils
from camel.app.tools.pipelines.mycobacterium.bamaddcustomtag import BAMAddCustomTag
from camel.scripts.mycobacteriumpipeline import CONFIG_DATA
from camel.scripts.mycobacteriumpipeline.mainmycobacteriumpipeline import main
from camel.tests import longRunningTest

CUSTOM_ANALYSES = basepipeutils.get_custom_analyses(CONFIG_DATA)


class TestMycobacteriumPipeline(CamelTestSuite):
    """
    Tests for the Mycobacterium pipeline.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir()
    input_fastq_pe = [
        test_file_dir / 'pipelines' / 'Myco-DRR041783-ds_1.fastq.gz',
        test_file_dir / 'pipelines' / 'Myco-DRR041783-ds_2.fastq.gz',
    ]
    input_fasta = test_file_dir / 'pipelines' / 'Myco-DRR041783-ds.fasta'
    input_fasta_csbrd = test_file_dir / 'pipelines' / 'Myco-S25BD09545.fasta'
    input_fastq_se = test_file_dir / 'pipelines' / 'Myco-SRR8948399_ont-ds.fastq.gz'
    input_vcf = test_file_dir / 'pipelines' / 'variants-Myco-DRR041783-ds-all.vcf'

    @longRunningTest()
    def test_blast(self) -> None:
        """
        Tests the Mycobacterium pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        path_bam_out = self.running_dir / 'out' / 'mapping.bam'
        result = cliutils.invoke(
            main,
            [
                '--fastq-pe',
                str(TestMycobacteriumPipeline.input_fastq_pe[0]),
                str(TestMycobacteriumPipeline.input_fastq_pe[1]),
                '--input-type',
                'illumina',
                '--output-html',
                str(path_report_out),
                '--output-dir',
                str(path_report_out.parent),
                '--output-tsv',
                str(path_summary_out),
                '--output-bam',
                str(path_bam_out),
                '--working-dir',
                str(self.running_dir),
                '--analyses',
                ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
                '--threads',
                '4',
            ],
        )
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)
        self.assertGreater(path_bam_out.stat().st_size, 0)

    @longRunningTest()
    def test_kma(self) -> None:
        """
        Tests the Mycobacterium pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(
            main,
            [
                '--fastq-pe',
                str(TestMycobacteriumPipeline.input_fastq_pe[0]),
                str(TestMycobacteriumPipeline.input_fastq_pe[1]),
                '--input-type',
                'illumina',
                '--output-html',
                str(path_report_out),
                '--output-dir',
                str(path_report_out.parent),
                '--output-tsv',
                str(path_summary_out),
                '--working-dir',
                str(self.running_dir),
                '--typing-method',
                'kma',
                '--gene-detection-method',
                'kma',
                '--analyses',
                ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
                '--threads',
                '4',
            ],
        )
        self.assertEqual(result.exit_code, 0)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_fasta_csbrd(self) -> None:
        """
        Tests the Mycobacterium pipeline using FASTA as input with only the rdcsb assay as a short running test.
        The rdcsb reporter relies on an external apt package and the error of the missing package was initially not
        caught upon migration to noble.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(
            main,
            [
                '--fasta',
                str(TestMycobacteriumPipeline.input_fasta_csbrd),
                '--input-type',
                'fasta',
                '--output-html',
                str(path_report_out),
                '--output-dir',
                str(path_report_out.parent),
                '--output-tsv',
                str(path_summary_out),
                '--working-dir',
                str(self.running_dir),
                '--analyses',
                'csb_rd',
                '--threads',
                '4',
            ],
        )
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_fasta(self) -> None:
        """
        Tests the Mycobacterium pipeline using FASTA as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(
            main,
            [
                '--fasta',
                str(TestMycobacteriumPipeline.input_fasta),
                '--input-type',
                'fasta',
                '--output-html',
                str(path_report_out),
                '--output-dir',
                str(path_report_out.parent),
                '--output-tsv',
                str(path_summary_out),
                '--working-dir',
                str(self.running_dir),
                '--analyses',
                ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
                '--threads',
                '4',
            ],
        )
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_fasta_with_vcf(self) -> None:
        """
        Tests the Mycobacterium pipeline using FASTA with VCF as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(
            main,
            [
                '--fasta',
                str(TestMycobacteriumPipeline.input_fasta),
                '--vcf-unfiltered',
                str(TestMycobacteriumPipeline.input_vcf),
                '--input-type',
                'fasta_with_vcf',
                '--output-html',
                str(path_report_out),
                '--output-dir',
                str(path_report_out.parent),
                '--output-tsv',
                str(path_summary_out),
                '--working-dir',
                str(self.running_dir),
                '--analyses',
                ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
                '--threads',
                '4',
            ],
        )
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_ont(self) -> None:
        """
        Tests the Mycobacterium pipeline using ONT as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(
            main,
            [
                '--fastq-se',
                str(TestMycobacteriumPipeline.input_fastq_se),
                '--input-type',
                'ont',
                '--output-html',
                str(path_report_out),
                '--output-dir',
                str(path_report_out.parent),
                '--output-tsv',
                str(path_summary_out),
                '--working-dir',
                str(self.running_dir),
                '--analyses',
                ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
                '--threads',
                '4',
            ],
        )
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_kma_ont(self) -> None:
        """
        Tests the Mycobacterium pipeline with all assays except for cgMLST, ONT input and kma as detection method
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(
            main,
            [
                '--fastq-se',
                str(TestMycobacteriumPipeline.input_fastq_se),
                '--input-type',
                'ont',
                '--output-html',
                str(path_report_out),
                '--output-dir',
                str(path_report_out.parent),
                '--output-tsv',
                str(path_summary_out),
                '--working-dir',
                str(self.running_dir),
                '--typing-method',
                'kma',
                '--gene-detection-method',
                'kma',
                '--analyses',
                ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
                '--threads',
                '4',
            ],
        )
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_add_custom_tag(self) -> None:
        """
        Tests the tool that adds a custom tag to the BAM output (used for using the BAM output in PACU).
        :return: None
        """
        path_test = TestMycobacteriumPipeline.test_file_dir / 'components' / 'toy.bam'
        add_tag = BAMAddCustomTag()
        add_tag.add_input_files({'BAM': [ToolIOFile(path_test)]})
        add_tag.update_parameters(output='custom_tag.bam', name='CT', value='my_value')
        add_tag.run(self.running_dir)
        self.verify_output_files(add_tag, 'BAM')


if __name__ == '__main__':
    unittest.main()
