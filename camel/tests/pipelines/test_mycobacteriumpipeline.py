import unittest
from pathlib import Path

import yaml


from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.pipelines.mycobacterium.bamaddcustomtag import BAMAddCustomTag
from camel.scripts.mycobacteriumpipeline import CONFIG_DATA
from camel.scripts.mycobacteriumpipeline.mainmycobacteriumpipeline import MainMycobacteriumPipeline
from camel.tests import longRunningTest


class TestMycobacteriumPipeline(CamelTestSuite):
    """
    Tests for the Mycobacterium pipeline.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir()
    input_fastq_pe = [
        test_file_dir / 'pipelines' / 'Myco-DRR041783-ds_1.fastq.gz',
        test_file_dir / 'pipelines' / 'Myco-DRR041783-ds_2.fastq.gz'
    ]
    input_fasta = test_file_dir / 'pipelines' / 'Myco-DRR041783-ds.fasta'
    input_fastq_se = test_file_dir / 'pipelines' / 'Myco-SRR8948399_ont-ds.fastq.gz'
    input_vcf = test_file_dir / 'pipelines' / 'variants-Myco-DRR041783-ds-all.vcf'

    def test_mycobacterium_pipeline_typing_db(self) -> None:
        """
        Checks if the databases for the sequence typing are available.
        :return: None
        """
        from camel.app.tools.pipelines.sequence_typing.locussetmanager import LocusSetManager
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
    def test_mycobacterium_pipeline_blast(self) -> None:
        """
        Tests the Mycobacterium pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        path_bam_out = self.running_dir / 'out' / 'mapping.bam'
        args = [
            '--fastq-pe', str(TestMycobacteriumPipeline.input_fastq_pe[0]),
            str(TestMycobacteriumPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--output-bam', str(path_bam_out),
            '--working-dir', str(self.running_dir)
        ] + [f"--{a.replace('_', '-')}" for a in MainMycobacteriumPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainMycobacteriumPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)
        self.assertGreater(path_bam_out.stat().st_size, 0)

    @longRunningTest()
    def test_mycobacterium_pipeline_kma(self) -> None:
        """
        Tests the Mycobacterium pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestMycobacteriumPipeline.input_fastq_pe[0]),
            str(TestMycobacteriumPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma'
        ] + [f"--{a.replace('_', '-')}" for a in MainMycobacteriumPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainMycobacteriumPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_mycobacterium_pipeline_fasta(self) -> None:
        """
        Tests the Mycobacterium pipeline using FASTA as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fasta', str(TestMycobacteriumPipeline.input_fasta),
                   '--input-type', 'fasta',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir)
               ] + [f"--{a.replace('_', '-')}" for a in MainMycobacteriumPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainMycobacteriumPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_mycobacterium_pipeline_fasta_with_vcf(self) -> None:
        """
        Tests the Mycobacterium pipeline using FASTA with VCF as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fasta', str(TestMycobacteriumPipeline.input_fasta),
                   '--vcf-unfiltered', str(TestMycobacteriumPipeline.input_vcf),
                   '--input-type', 'fasta_with_vcf',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir)
               ] + [f"--{a.replace('_', '-')}" for a in MainMycobacteriumPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainMycobacteriumPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_mycobacterium_pipeline_ont(self) -> None:
        """
        Tests the Mycobacterium pipeline using ONT as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fastq-se', str(TestMycobacteriumPipeline.input_fastq_se),
                   '--input-type', 'ont',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir)
               ] + [f"--{a.replace('_', '-')}" for a in MainMycobacteriumPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainMycobacteriumPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_mycobacterium_pipeline_kma_ont(self) -> None:
        """
        Tests the Mycobacterium pipeline with all assays except for cgMLST, ONT input and kma as detection method
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fastq-se', str(TestMycobacteriumPipeline.input_fastq_se),
                   '--input-type', 'ont',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir),
                   '--detection-method', 'kma'
               ] + [f"--{a.replace('_', '-')}" for a in MainMycobacteriumPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainMycobacteriumPipeline(args)
        main.run()
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
