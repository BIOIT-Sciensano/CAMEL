import unittest
from pathlib import Path

from camelcore.app.io.tooliodirectory import ToolIODirectory
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.config import config
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.kraken.kraken2 import Kraken2
from camel.tests import requires_dependency_service


class TestKraken2(CamelTestSuite):
    """
    Tests the Kraken2 tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('kraken')
    input_pe_reads = [test_file_dir / 'lm1_1.fastq', test_file_dir / 'lm1_2.fastq']
    input_fasta_reads = test_file_dir / 'reads_test.fasta'
    input_db = Path(config.dir_db, 'kraken2_microbial', 'latest')

    @requires_dependency_service('lmod')
    def test_kraken2_paired(self) -> None:
        """
        Tests kraken2 with paired-end input.
        :return: None
        """
        kraken2 = Kraken2()
        kraken2.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in self.input_pe_reads],
            'DB': [ToolIODirectory(TestKraken2.input_db)]
        })
        kraken2.run(self.running_dir)
        self.verify_output_files(kraken2, 'TSV')
        self.verify_output_files(kraken2, 'TSV_report')

    @requires_dependency_service('lmod')
    def test_kraken2_fasta(self) -> None:
        """
        Tests kraken2 with fasta input.
        :return: None
        """
        kraken2 = Kraken2()
        kraken2.add_input_files({
            'FASTA': [ToolIOFile(self.input_fasta_reads) ],
            'DB': [ToolIODirectory(TestKraken2.input_db)]
        })
        kraken2.run(self.running_dir)
        self.verify_output_files(kraken2, 'TSV')
        self.verify_output_files(kraken2, 'TSV_report')


if __name__ == '__main__':
    unittest.main()
