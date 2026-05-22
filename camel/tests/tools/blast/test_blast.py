import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.config import config
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.loggers import initialize_logging
from camel.app.tools.blast.blastformatter import BlastFormatter
from camel.app.tools.blast.blastn import Blastn
from camel.app.tools.blast.blastx import Blastx


class TestBlast(CamelTestSuite):
    """
    Tests for BLAST+.
    """

    def test_blastn(self) -> None:
        """
        Tests the blastn tool with default settings.
        :return: None
        """
        blastn = Blastn()
        blastn.add_input_files({
            'FASTA': [ToolIOFile(Path(config.dir_testdata, 'blast', 'plasmid_contig.fasta'))],
            'FASTA_Subject': [ToolIOFile(Path(config.dir_testdata, 'blast', 'amr_genes.fasta'))],
        })
        blastn.update_parameters(task='megablast', output_format='6', threads=1)
        blastn.run(self.running_dir)
        self.verify_output_files(blastn, 'TSV', 1)

    def test_blastn_asn_output(self) -> None:
        """
        Tests the blastn tool with ASN output.
        :return: None
        """
        blastn = Blastn()
        blastn.add_input_files({
            'FASTA': [ToolIOFile(Path(config.dir_testdata, 'blast', 'plasmid_contig.fasta'))],
            'FASTA_Subject': [ToolIOFile(Path(config.dir_testdata, 'blast', 'amr_genes.fasta'))],
        })
        blastn.update_parameters(task='megablast', output_format='11', output_filename='blast.asn', threads=1)
        blastn.run(self.running_dir)
        self.verify_output_files(blastn, 'ASN', 1)

    def test_blast_formatter(self) -> None:
        """
        Tests the blast formatter tool with ASN output from blastn.
        :return: None
        """
        # Run blastn
        blastn = Blastn()
        blastn.add_input_files({
            'FASTA': [ToolIOFile(Path(config.dir_testdata, 'blast', 'plasmid_contig.fasta'))],
            'FASTA_Subject': [ToolIOFile(Path(config.dir_testdata, 'blast', 'amr_genes.fasta'))],
        })
        blastn.update_parameters(task='megablast', output_format='11', output_filename='blast.asn', threads=1)
        blastn.run(self.running_dir)

        # Run blast formatter
        blast_formatter = BlastFormatter()
        blast_formatter.add_input_files({'ASN': blastn.tool_outputs['ASN']})
        blast_formatter.update_parameters(output_format='6')
        blast_formatter.run(self.running_dir)
        self.verify_output_files(blast_formatter, 'TSV', 1)

    def test_blastx(self) -> None:
        """
        Tests the blastx tool with default settings.
        :return: None
        """
        blastn = Blastx()
        blastn.add_input_files({
            'FASTA': [ToolIOFile(Path(config.dir_testdata, 'blast', 'plasmid_contig.fasta'))],
            'FASTA_Subject': [ToolIOFile(Path(config.dir_testdata, 'blast', 'PorA_VR1.fasta'))],
        })
        blastn.update_parameters(output_format='6', threads=1)
        blastn.run(self.running_dir)
        self.verify_output_files(blastn, 'TSV', 1)


if __name__ == '__main__':
    initialize_logging()
    unittest.main()
