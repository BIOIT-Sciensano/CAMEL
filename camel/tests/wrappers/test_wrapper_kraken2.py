import unittest
from pathlib import Path

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.config import config
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.scriptutils.basepipe.fastqinput import FastqInput
from camel.app.wrappers.kraken2wrapper import Kraken2Wrapper


class TestWrapperKraken2(CamelTestSuite):
    """
    Tests the Kraken2 contamination detection wrapper.
    """
    # Input files (Gene detection)
    test_file_dir = CamelTestSuite.get_test_file_dir('workflows')

    # Input files
    fastq_pe_contamination = [
        test_file_dir / 'kraken2' / 'reads_illumina_1.fastq',
        test_file_dir / 'kraken2' / 'reads_illumina_2.fastq']
    fastq_pe = [
        test_file_dir / 'kraken2' / 'neisseria_10k_1.fastq.gz',
        test_file_dir / 'kraken2' / 'neisseria_10k_2.fastq.gz']

    fastq_se_nanopore = test_file_dir / 'kraken2' / 'reads_nanopore_1.fastq.gz'
    fastq_se_nanopore_contamination = test_file_dir / 'kraken2' / 'reads_nanopore_1.fastq.gz'

    fasta = test_file_dir / 'contigs.fasta'
    path_db = Path(config.dir_db / 'kraken2_microbial' / 'latest')

    def test_kraken2_illumina_paired_end(self) -> None:
        """
        Tests the KRAKEN2 wrapper on Illumina PE data.
        :return: None
        """
        wrapper = Kraken2Wrapper(self.running_dir)
        fq_in = FastqInput('PE', pe=[ToolIOFile(fq) for fq in TestWrapperKraken2.fastq_pe])
        expected_species = 'Neisseria meningitidis'
        wrapper.run_fastq('test_sample', fq_in, 'illumina', expected_species, TestWrapperKraken2.path_db)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)
        self.assertGreater(wrapper.output.tsv_summary.stat().st_size, 0)
        CamelTestSuite.export_report_section(wrapper.output.report_section, self.running_dir / 'report')

    def test_kraken2_nanopore_single_end(self) -> None:
        """
        Tests the KRAKEN2 wrapper on Nanopore SE data.
        :return: None
        """
        wrapper = Kraken2Wrapper(self.running_dir)
        expected_species = 'Influenza A virus'
        fq_in = FastqInput('SE', se=[ToolIOFile(TestWrapperKraken2.fastq_se_nanopore)], is_pe=False)
        wrapper.run_fastq('test_sample', fq_in, 'ont', expected_species, db=TestWrapperKraken2.path_db)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)
        self.assertGreater(wrapper.output.tsv_summary.stat().st_size, 0)
        CamelTestSuite.export_report_section(wrapper.output.report_section, self.running_dir / 'report')

    def test_kraken2_illumina_paired_end_contaminated(self) -> None:
        """
        Tests the KRAKEN2 wrapper on Illumina PE data, with reads other than the expected species.
        :return: None
        """
        wrapper = Kraken2Wrapper(self.running_dir)
        fq_in = FastqInput('PE', pe=[ToolIOFile(fq) for fq in TestWrapperKraken2.fastq_pe_contamination])
        expected_species = 'Escherichia coli'
        wrapper.run_fastq('test_sample', fq_in, 'illumina', expected_species, db=TestWrapperKraken2.path_db)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)
        self.assertGreater(wrapper.output.tsv_summary.stat().st_size, 0)
        list_contaminants_to_test = \
            wrapper.output.informs['contaminants_warn'] + wrapper.output.informs['contaminants_fail']
        self.assertGreater(len(list_contaminants_to_test), 0)
        CamelTestSuite.export_report_section(wrapper.output.report_section, self.running_dir / 'report')

    def test_kraken2_nanopore_single_end_contaminated(self) -> None:
        """
        Tests the KRAKEN2 wrapper on Nanopore SE data, with reads other than the expected species.
        :return: None
        """
        wrapper = Kraken2Wrapper(self.running_dir)
        expected_species = 'Influenza A virus'
        fq_in = FastqInput('SE', se=[ToolIOFile(TestWrapperKraken2.fastq_se_nanopore_contamination)], is_pe=False)
        wrapper.run_fastq('test_sample', fq_in, 'ont', expected_species, db=TestWrapperKraken2.path_db)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)
        self.assertGreater(wrapper.output.tsv_summary.stat().st_size, 0)
        CamelTestSuite.export_report_section(wrapper.output.report_section, self.running_dir / 'report')

    def test_kraken2_illumina_paired_end_genus(self) -> None:
        """
        Tests the KRAKEN2 wrapper on Illumina PE data.
        :return: None
        """
        wrapper = Kraken2Wrapper(self.running_dir)
        expected_genus = 'Neisseria'
        fq_in = FastqInput('PE', pe=[ToolIOFile(fq) for fq in TestWrapperKraken2.fastq_pe])
        wrapper.run_fastq(
            'test_sample', fq_in, 'illumina', expected_genus, db=TestWrapperKraken2.path_db, level_of_depth='G')
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)
        self.assertGreater(wrapper.output.tsv_summary.stat().st_size, 0)
        CamelTestSuite.export_report_section(wrapper.output.report_section, self.running_dir / 'report')

    def test_kraken2_illumina_paired_end_genus_contamination(self) -> None:
        """
        Tests the KRAKEN2 wrapper on Illumina PE data.
        :return: None
        """
        wrapper = Kraken2Wrapper(self.running_dir)
        expected_genus = 'Escherichia'
        fq_in = FastqInput('PE', pe=[ToolIOFile(fq) for fq in TestWrapperKraken2.fastq_pe_contamination])
        wrapper.run_fastq(
            'test_sample', fq_in, 'illumina', expected_genus, db=TestWrapperKraken2.path_db, level_of_depth='G')
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)
        self.assertGreater(wrapper.output.tsv_summary.stat().st_size, 0)
        list_contaminants_to_test = \
            wrapper.output.informs['contaminants_warn'] + wrapper.output.informs['contaminants_fail']
        self.assertGreater(len(list_contaminants_to_test), 0)
        CamelTestSuite.export_report_section(wrapper.output.report_section, self.running_dir / 'report')

    def test_kraken2_fasta(self) -> None:
        """
        Tests the KRAKEN2 wrapper on FASTA data.
        :return: None
        """
        wrapper = Kraken2Wrapper(self.running_dir)
        expected_species = 'Escherichia coli'
        wrapper.run_fasta(
            'test_sample', TestWrapperKraken2.fasta, expected_species, db=TestWrapperKraken2.path_db)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)
        self.assertGreater(wrapper.output.tsv_summary.stat().st_size, 0)
        CamelTestSuite.export_report_section(wrapper.output.report_section, self.running_dir / 'report')

    def test_kraken2_fasta_contamination(self) -> None:
        """
        Tests the KRAKEN2 wrapper on FASTA data not from the expected species.
        :return: None
        """
        wrapper = Kraken2Wrapper(self.running_dir)
        expected_species = 'Neisseria meningitidis'
        wrapper.run_fasta(
            'test_sample', TestWrapperKraken2.fasta, expected_species, db=TestWrapperKraken2.path_db)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)
        self.assertGreater(wrapper.output.tsv_summary.stat().st_size, 0)
        list_contaminants_to_test = \
            wrapper.output.informs['contaminants_warn'] + wrapper.output.informs['contaminants_fail']
        self.assertGreater(len(list_contaminants_to_test), 0)
        CamelTestSuite.export_report_section(wrapper.output.report_section, self.running_dir / 'report')


if __name__ == '__main__':
    unittest.main()
