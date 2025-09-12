import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.pipelines.salmonella.seqsero2 import SeqSero2
from camel.app.tools.pipelines.salmonella.seqsero2reporter import SeqSero2Reporter


class TestSeqsero2(CamelTestSuite):
    """
    Tests the Seqsero2 tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('salmonella')
    input_pe_reads = [test_file_dir / "SRR493330_1.fastq.gz", test_file_dir / "SRR493330_2.fastq.gz"]
    input_fasta_file = test_file_dir / 'assembly_filtered.fasta'
    input_fastq_se = test_file_dir / 'Salmonella_S23BD05337-RBK_ont-ds-ds.fastq.gz'
    db_path = Path(CamelTestSuite.camel.config['db_root']) / 'pipelines/salmonella/seqsero2/1.2.1/seqsero2_db'

    def test_seqsero2_kmer(self) -> None:
        """
        Tests basic seqsero2 run in kmer mode.
        :return: None
        """
        seqsero2_tool = SeqSero2()
        seqsero2_tool.add_input_files({
            'FASTA': [ToolIOFile(Path(self.input_fasta_file))],
            'DIR': [ToolIODirectory(self.db_path)]
        })
        seqsero2_tool.update_parameters(mode='kmer')
        seqsero2_tool.run(self.running_dir)
        self.verify_output_files(seqsero2_tool, 'TXT')
        self.assertIn('db_path', seqsero2_tool.informs)

    def test_seqsero2_allele(self) -> None:
        """
        Tests basic seqsero2 run in allele mode.
        :return: None
        """
        seqsero2_tool = SeqSero2()
        seqsero2_tool.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in self.input_pe_reads],
            'DIR': [ToolIODirectory(self.db_path)]
        })
        seqsero2_tool.update_parameters(mode='allele')
        seqsero2_tool.run(self.running_dir)
        self.verify_output_files(seqsero2_tool, 'TXT')
        self.assertIn('db_path', seqsero2_tool.informs)

    def test_seqsero2_kmerread(self) -> None:
        """
        Tests basic seqsero2 run in kmerread mode.
        :return: None
        """
        seqsero2_tool = SeqSero2()
        seqsero2_tool.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in self.input_pe_reads],
            'DIR': [ToolIODirectory(self.db_path)]
        })
        seqsero2_tool.update_parameters(mode='kmerread')
        seqsero2_tool.run(self.running_dir)
        self.verify_output_files(seqsero2_tool, 'TXT')
        self.assertIn('db_path', seqsero2_tool.informs)

    def test_seqsero2_ont(self) -> None:
        """
        Tests basic seqsero2 run in kmerread mode with ont data input
        :return: None
        """
        seqsero2_tool = SeqSero2()
        seqsero2_tool.add_input_files({
            'FASTQ_ONT': [ToolIOFile(self.input_fastq_se)],
            'DIR': [ToolIODirectory(self.db_path)]
        })
        seqsero2_tool.update_parameters(mode='kmerread')
        seqsero2_tool.run(self.running_dir)
        self.verify_output_files(seqsero2_tool, 'TXT')
        self.assertIn('db_path', seqsero2_tool.informs)

    def test_seqsero2_reporter(self) -> None:
        """
        Tests seqsero2 reporter.
        :return: None
        """
        seqsero2_tool = SeqSero2()
        seqsero2_tool.add_input_files({
            'FASTA': [ToolIOFile(Path(self.input_fasta_file))],
            'DIR': [ToolIODirectory(self.db_path)]
        })
        seqsero2_tool.update_parameters(mode='kmer')
        seqsero2_tool.run(self.running_dir)
        self.verify_output_files(seqsero2_tool, 'TXT')
        self.assertIn('db_path', seqsero2_tool.informs)

        seqsero2_reporter = SeqSero2Reporter()
        seqsero2_reporter.add_input_files({'TXT_seqsero2_kmer': seqsero2_tool.tool_outputs['TXT'],
                                           'DIR_seqsero2': [ToolIODirectory(self.db_path)]})
        seqsero2_reporter.add_input_informs({'serotyping_seqsero2': seqsero2_tool.informs})
        seqsero2_reporter.run(self.running_dir)
        output_section = seqsero2_reporter.tool_outputs['VAL_HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
