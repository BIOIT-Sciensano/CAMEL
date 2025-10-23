import unittest
from pathlib import Path

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.pipelines.salmonella.spifinder import SPIFinder
from camel.app.tools.pipelines.salmonella.spifinderreporter import SPIFinderReporter


class TestSPIFinder(CamelTestSuite):
    """
    Tests the SPIFinder tool and its reporter.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('salmonella')
    input_pe_reads = [test_file_dir / "SRR493330_1.fastq.gz", test_file_dir / "SRR493330_2.fastq.gz"]
    input_fasta_file = test_file_dir / 'assembly_filtered.fasta'
    db_path = Path('/db/pipelines/salmonella/spifinder/genomicepidemiology-spifinder_db-db102668b704')

    def test_spifinder_fastq(self) -> None:
        """
        Tests basic spifinder run on fastq input data.
        :return: None
        """
        spifinder = SPIFinder()
        spifinder.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in self.input_pe_reads],
            'DIR': [ToolIODirectory(self.db_path)]
        })
        spifinder.run(self.running_dir)
        self.verify_output_files(spifinder, 'JSON')

    def test_spifinder_fasta(self) -> None:
        """
        Tests basic spifinder run on fasta input data.
        :return: None
        """
        spifinder = SPIFinder()
        spifinder.add_input_files({
            'FASTA': [ToolIOFile(Path(TestSPIFinder.input_fasta_file))],
            'DIR': [ToolIODirectory(self.db_path)]
        })
        spifinder.run(self.running_dir)
        self.verify_output_files(spifinder, 'JSON')

    def test_spifinder_reporter(self) -> None:
        """
        Tests the spifinder reporter tool.
        :return: None
        """
        spifinder_fastq = SPIFinder()
        spifinder_fastq.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in self.input_pe_reads],
            'DIR': [ToolIODirectory(self.db_path)]
        })
        spifinder_fastq.run(self.running_dir)
        self.verify_output_files(spifinder_fastq, 'JSON')

        spifinder_fasta = SPIFinder()
        spifinder_fasta.add_input_files({
            'FASTA': [ToolIOFile(Path(TestSPIFinder.input_fasta_file))],
            'DIR': [ToolIODirectory(self.db_path)]
        })
        spifinder_fasta.run(self.running_dir)
        self.verify_output_files(spifinder_fasta, 'JSON')

        dummy_tsv_spifinder_documentation = Path('./dummy.tsv')
        dummy_tsv_spifinder_documentation.touch()
        spifinderreporter = SPIFinderReporter()
        spifinderreporter.add_input_files({'TSV_documentation': [ToolIOFile(dummy_tsv_spifinder_documentation)],
                                           'JSON_FASTQ': spifinder_fastq.tool_outputs['JSON'],
                                           'JSON_FASTA': spifinder_fasta.tool_outputs['JSON']})
        spifinderreporter.add_input_informs({'spifinder_fastq': spifinder_fastq.informs,
                                             'spifinder_fasta': spifinder_fasta.informs})
        spifinderreporter.run(self.running_dir)
        self.assertGreater(len(spifinderreporter.tool_outputs['VAL_HTML'][0].value.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
