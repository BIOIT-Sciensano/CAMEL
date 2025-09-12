import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.resfinder.resfinder import ResFinder
from camel.app.tools.resfinder.resfinderreporter import ResFinderReporter


class TestResFinder(CamelTestSuite):
    """
    Initializes this testing tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('resfinder')
    FILE_FASTA_1 = ToolIOFile(test_file_dir / 'ref_ecoli.fasta')
    FILE_FASTA_2 = ToolIOFile(test_file_dir / 'salmonella_lt2_ref.fasta')
    FILE_FASTA_3 = ToolIOFile(test_file_dir / 'assembly-VAR305.fasta')
    FILE_FASTA_4 = ToolIOFile(test_file_dir / 'enterococcus_unkown_amr_muts.fasta')
    FILE_FASTQ_1 = ToolIOFile(test_file_dir / 'reads_illumina_1.fastq')
    FILE_FASTQ_2 = ToolIOFile(test_file_dir / 'reads_illumina_2.fastq')
    FILE_FASTA_ENTERO = ToolIOFile(test_file_dir / 'assembly_entero.fasta')
    DB_RESFINDER = Path(CamelTestSuite.camel.config['db_root'], 'resfinder4')

    def test_resfinder_fasta(self) -> None:
        """
        Actually testing ResFinder with contigs file.
        :return: None
        """
        resfinder = ResFinder()
        resfinder.add_input_files({'FASTA': [TestResFinder.FILE_FASTA_1], 'DIR': [ToolIODirectory(self.DB_RESFINDER)]})
        resfinder.update_parameters(output_path=self.running_dir, min_cov=0.6, threshold=0.8, acquired=True)
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV_genes')
        self.verify_output_files(resfinder, 'TSV_pheno_general')

    def test_resfinder_fastq(self) -> None:
        """
        Testing resfinder with paired-end fastq reads.
        :return: None
        """
        resfinder = ResFinder()
        resfinder.add_input_files({
            'FASTQ_PE': [TestResFinder.FILE_FASTQ_1, TestResFinder.FILE_FASTQ_2],
            'DIR': [ToolIODirectory(self.DB_RESFINDER)]})
        resfinder.update_parameters(output_path=self.running_dir, min_cov=0.6, threshold=0.8, acquired=True)
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV_genes')
        self.verify_output_files(resfinder, 'TSV_pheno_general')

    def test_resfinder_pointfinder_fasta(self) -> None:
        """
        Testing resfinder with pointfinder mode and fasta file.
        :return: None
        """
        resfinder = ResFinder()
        resfinder.add_input_files({
            'FASTA': [TestResFinder.FILE_FASTA_1],
            'DIR': [ToolIODirectory(self.DB_RESFINDER)]})
        resfinder.update_parameters(
            output_path=self.running_dir, min_cov=0.6, threshold=0.8, point=True, species='"escherichia coli"')
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV_point')
        self.verify_output_files(resfinder, 'TSV_pheno_general')
        self.verify_output_files(resfinder, 'TSV_pheno_species')

    def test_resfinder_pointfinder_fastq(self) -> None:
        """
        Testing resfinder with pointfinder mode and fastq files.
        :return: None
        """
        resfinder = ResFinder()
        resfinder.add_input_files({
            'FASTQ_PE': [TestResFinder.FILE_FASTQ_1, TestResFinder.FILE_FASTQ_2],
            'DIR': [ToolIODirectory(self.DB_RESFINDER)]})
        resfinder.update_parameters(
            output_path=self.running_dir, min_cov=0.6, threshold=0.8, point=True, species='"escherichia coli"')
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV_point')
        self.verify_output_files(resfinder, 'TSV_pheno_general')
        self.verify_output_files(resfinder, 'TSV_pheno_species')

    def test_resfinder_fasta_enterococcus(self) -> None:
        """
        Testing resfinder with an Enterococcus input FASTA file.
        :return: None
        """
        resfinder = ResFinder()
        resfinder.add_input_files({
            'FASTA': [TestResFinder.FILE_FASTA_ENTERO],
            'DIR': [ToolIODirectory(self.DB_RESFINDER)]
        })
        resfinder.update_parameters(
            output_path=self.running_dir, min_cov=0.9, threshold=0.8, point=True, acquired=True,
            species='"enterococcus faecalis"')
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV_point')
        self.verify_output_files(resfinder, 'TSV_genes')
        self.verify_output_files(resfinder, 'TSV_pheno_general')
        self.verify_output_files(resfinder, 'TSV_pheno_species')

        reporter = ResFinderReporter()
        reporter.add_input_files({
            'TSV_point': resfinder.tool_outputs['TSV_point'],
            'TSV_pheno_general': resfinder.tool_outputs['TSV_pheno_general'],
            'TSV_pheno_species': resfinder.tool_outputs['TSV_pheno_species'],
            'TSV_genes': resfinder.tool_outputs['TSV_genes']})
        reporter.add_input_informs({'resfinder': resfinder.informs})
        reporter.run(self.running_dir)
        output_section = reporter.tool_outputs['VAL_HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)

    def test_resfinder_fasta_enterococcus_no_species(self) -> None:
        """
        Testing resfinder with an Enterococcus input FASTA file.
        :return: None
        """
        # Run ResFinder
        resfinder = ResFinder()
        resfinder.add_input_files({
            'FASTA': [TestResFinder.FILE_FASTA_ENTERO],
            'DIR': [ToolIODirectory(self.DB_RESFINDER)]
        })
        resfinder.update_parameters(
            output_path=self.running_dir, min_cov=0.9, threshold=0.8, point=False, acquired=True)
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV_genes')
        self.verify_output_files(resfinder, 'TSV_pheno_general')

        # Run the reporter
        reporter = ResFinderReporter()
        reporter.add_input_files({
            'TSV_pheno_general': resfinder.tool_outputs['TSV_pheno_general'],
            'TSV_genes': resfinder.tool_outputs['TSV_genes']})
        reporter.add_input_informs({'resfinder': resfinder.informs})
        reporter.run(self.running_dir)
        output_section = reporter.tool_outputs['VAL_HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
