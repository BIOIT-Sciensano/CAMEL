import logging
import unittest
from pathlib import Path

from camelcore.app.io.tooliodirectory import ToolIODirectory
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.io.tooliovalue import ToolIOValue

from camel.app.config import config
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mykrobe.mykrobe import Mykrobe
from camel.app.tools.mykrobe.mykrobereporter import MykrobeReporter


class TestMykrobe(CamelTestSuite):
    """
    Tests the Mykrobe tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('salmonella')
    testdata_dir = Path('/testdata/camel/pipelines')
    db_dir = Path(config.dir_db, 'mykrobe/latest')

    typhi_pe_reads = [
        test_file_dir / "SRR493330_1.fastq.gz",
        test_file_dir / "SRR493330_2.fastq.gz",
    ]
    typhi_ont = [testdata_dir / 'Salmonella_ERR11177482.fastq.gz']

    shigella_pe_reads = [
        testdata_dir / 'Shigella-S17BD07654_1.fastq.gz',
        testdata_dir / 'Shigella-S17BD07654_2.fastq.gz',
    ]
    shigella_fasta = [testdata_dir / 'Shigella-S17BD07654.fasta']

    staph_pe_reads = [
        testdata_dir / 'Saureus-SRR10393587-ds_1.fastq.gz',
        testdata_dir / 'Saureus-SRR10393587-ds_2.fastq.gz',
    ]

    myco_pe_reads = [
        testdata_dir / 'Myco-DRR041783-ds_1.fastq.gz',
        testdata_dir / 'Myco-DRR041783-ds_2.fastq.gz',
    ]
    # Output file
    output_csv = test_file_dir / 'mykrobe.csv'

    def test_typhi_mykrobe(self) -> None:
        """
        Tests Mykrobe on Salmonella Typhi FASTQ_PE reads.
        :return: None
        """
        tool = Mykrobe()
        tool.add_input_files(
            {
                'FASTQ_PE': [ToolIOFile(x) for x in self.typhi_pe_reads],
                'DIR': [ToolIODirectory(self.db_dir)],
                'SPECIES': [ToolIOValue('typhi')],
            }
        )
        tool.run(self.running_dir)
        self.verify_output_files(tool, 'CSV')

    def test_typhi_ont_mykrobe(self) -> None:
        """
        Tests Mykrobe on Salmonella long-reads based FASTQ file.
        :return: None
        """
        tool = Mykrobe()
        tool.add_input_files(
            {
                'FASTQ_SE': [ToolIOFile(x) for x in self.typhi_ont],
                'DIR': [ToolIODirectory(self.db_dir)],
                'SPECIES': [ToolIOValue('typhi')],
            }
        )
        tool.run(self.running_dir)
        self.verify_output_files(tool, 'CSV')

    def test_shigella_pe_mykrobe(self) -> None:
        """
        Tests Mykrobe on Shigella FASTQ_PE reads.
        :return: None
        """
        tool = Mykrobe()
        tool.add_input_files(
            {
                'FASTQ_PE': [ToolIOFile(x) for x in self.shigella_pe_reads],
                'DIR': [ToolIODirectory(self.db_dir)],
                'SPECIES': [ToolIOValue('sonnei')],
            }
        )
        tool.run(self.running_dir)
        self.verify_output_files(tool, 'CSV')

    def test_shigella_fasta_mykrobe(self) -> None:
        """
        Tests Mykrobe on Shigella FASTA file.
        :return: None
        """
        tool = Mykrobe()
        tool.add_input_files(
            {
                'FASTA': [ToolIOFile(x) for x in self.shigella_fasta],
                'DIR': [ToolIODirectory(self.db_dir)],
                'SPECIES': [ToolIOValue('sonnei')],
            }
        )
        tool.run(self.running_dir)
        self.verify_output_files(tool, 'CSV')

    def test_staph_mykrobe(self) -> None:
        """
        Tests Mykrobe on Staphylococcus aureus FASTQ_PE reads.
        :return: None
        """
        tool = Mykrobe()
        tool.add_input_files(
            {
                'FASTQ_PE': [ToolIOFile(x) for x in self.staph_pe_reads],
                'DIR': [ToolIODirectory(self.db_dir)],
                'SPECIES': [ToolIOValue('staph')],
            }
        )
        tool.run(self.running_dir)
        self.verify_output_files(tool, 'CSV')

    def test_myco_mykrobe(self) -> None:
        """
        Tests Mykrobe on Mycobacterium tuberculosis FASTQ_PE reads.
        :return: None
        """
        tool = Mykrobe()
        tool.add_input_files(
            {
                'FASTQ_PE': [ToolIOFile(x) for x in self.myco_pe_reads],
                'DIR': [ToolIODirectory(self.db_dir)],
                'SPECIES': [ToolIOValue('tb')],
            }
        )
        tool.run(self.running_dir)
        self.verify_output_files(tool, 'CSV')

    def test_salmonella_reporter(self) -> None:
        """
        Tests Mykrobe reporter tool
        :return: None
        """
        # Run Mykrobe
        mykrobe = Mykrobe()
        mykrobe.add_input_files(
            {
                'FASTQ_PE': [ToolIOFile(x) for x in self.typhi_pe_reads],
                'DIR': [ToolIODirectory(self.db_dir)],
                'SPECIES': [ToolIOValue('typhi')],
            }
        )
        mykrobe.run(self.running_dir)
        self.verify_output_files(mykrobe, 'CSV')

        # Run the reporter
        mykrobereporter = MykrobeReporter()
        mykrobereporter.add_input_files({'CSV': mykrobe.tool_outputs['CSV']})
        mykrobereporter.add_input_informs({'mykrobe': mykrobe.informs})
        mykrobereporter.run(self.running_dir)
        self.assertGreater(
            len(mykrobereporter.tool_outputs['HTML'][0].value.to_html()), 0
        )

        # Save the report
        html_out = self.running_dir / 'report.html'
        with html_out.open('w') as handle:
            handle.write(mykrobereporter.tool_outputs['HTML'][0].value.to_html())
        logging.info(f'Output report created: {html_out}')

    def test_typhi_ont_reporter(self) -> None:
        """
        Tests Mykrobe reporter tool
        :return: None
        """
        # Run Mykrobe
        mykrobe = Mykrobe()
        mykrobe.add_input_files(
            {
                'FASTQ_SE': [ToolIOFile(x) for x in self.typhi_ont],
                'DIR': [ToolIODirectory(self.db_dir)],
                'SPECIES': [ToolIOValue('typhi')],
            }
        )
        mykrobe.run(self.running_dir)
        self.verify_output_files(mykrobe, 'CSV')

        # Run the reporter
        mykrobereporter = MykrobeReporter()
        mykrobereporter.add_input_files({'CSV': mykrobe.tool_outputs['CSV']})
        mykrobereporter.add_input_informs({'mykrobe': mykrobe.informs})
        mykrobereporter.run(self.running_dir)
        self.assertGreater(
            len(mykrobereporter.tool_outputs['HTML'][0].value.to_html()), 0
        )

        # Save the report
        html_out = self.running_dir / 'report.html'
        with html_out.open('w') as handle:
            handle.write(mykrobereporter.tool_outputs['HTML'][0].value.to_html())
        logging.info(f'Output report created: {html_out}')

    def test_shigella_reporter(self) -> None:
        """
        Tests Mykrobe reporter tool
        :return: None
        """
        # Run Mykrobe
        mykrobe = Mykrobe()
        mykrobe.add_input_files(
            {
                'FASTQ_PE': [ToolIOFile(x) for x in self.shigella_pe_reads],
                'DIR': [ToolIODirectory(self.db_dir)],
                'SPECIES': [ToolIOValue('sonnei')],
            }
        )
        mykrobe.run(self.running_dir)
        self.verify_output_files(mykrobe, 'CSV')

        # Run the reporter
        mykrobereporter = MykrobeReporter()
        mykrobereporter.add_input_files({'CSV': mykrobe.tool_outputs['CSV']})
        mykrobereporter.add_input_informs({'mykrobe': mykrobe.informs})
        mykrobereporter.run(self.running_dir)
        self.assertGreater(
            len(mykrobereporter.tool_outputs['HTML'][0].value.to_html()), 0
        )

        # Save the report
        html_out = self.running_dir / 'report.html'
        with html_out.open('w') as handle:
            handle.write(mykrobereporter.tool_outputs['HTML'][0].value.to_html())
        logging.info(f'Output report created: {html_out}')

    def test_shigella_fasta_reporter(self) -> None:
        """
        Tests Mykrobe reporter tool
        :return: None
        """
        # Run Mykrobe
        mykrobe = Mykrobe()
        mykrobe.add_input_files(
            {
                'FASTA': [ToolIOFile(x) for x in self.shigella_fasta],
                'DIR': [ToolIODirectory(self.db_dir)],
                'SPECIES': [ToolIOValue('sonnei')],
            }
        )
        mykrobe.run(self.running_dir)
        self.verify_output_files(mykrobe, 'CSV')

        # Run the reporter
        mykrobereporter = MykrobeReporter()
        mykrobereporter.add_input_files({'CSV': mykrobe.tool_outputs['CSV']})
        mykrobereporter.add_input_informs({'mykrobe': mykrobe.informs})
        mykrobereporter.run(self.running_dir)
        self.assertGreater(
            len(mykrobereporter.tool_outputs['HTML'][0].value.to_html()), 0
        )

        # Save the report
        html_out = self.running_dir / 'report.html'
        with html_out.open('w') as handle:
            handle.write(mykrobereporter.tool_outputs['HTML'][0].value.to_html())
        logging.info(f'Output report created: {html_out}')

    def test_staph_reporter(self) -> None:
        """
        Tests Mykrobe reporter tool
        :return: None
        """
        # Run Mykrobe
        mykrobe = Mykrobe()
        mykrobe.add_input_files(
            {
                'FASTQ_PE': [ToolIOFile(x) for x in self.staph_pe_reads],
                'DIR': [ToolIODirectory(self.db_dir)],
                'SPECIES': [ToolIOValue('staph')],
            }
        )
        mykrobe.run(self.running_dir)
        self.verify_output_files(mykrobe, 'CSV')

        # Run the reporter
        mykrobereporter = MykrobeReporter()
        mykrobereporter.add_input_files({'CSV': mykrobe.tool_outputs['CSV']})
        mykrobereporter.add_input_informs({'mykrobe': mykrobe.informs})
        mykrobereporter.run(self.running_dir)
        self.assertGreater(
            len(mykrobereporter.tool_outputs['HTML'][0].value.to_html()), 0
        )

        # Save the report
        html_out = self.running_dir / 'report.html'
        with html_out.open('w') as handle:
            handle.write(mykrobereporter.tool_outputs['HTML'][0].value.to_html())
        logging.info(f'Output report created: {html_out}')

    def test_myco_reporter(self) -> None:
        """
        Tests Mykrobe reporter tool
        :return: None
        """
        # Run Mykrobe
        mykrobe = Mykrobe()
        mykrobe.add_input_files(
            {
                'FASTQ_PE': [ToolIOFile(x) for x in self.myco_pe_reads],
                'DIR': [ToolIODirectory(self.db_dir)],
                'SPECIES': [ToolIOValue('tb')],
            }
        )
        mykrobe.run(self.running_dir)
        self.verify_output_files(mykrobe, 'CSV')

        # Run the reporter
        mykrobereporter = MykrobeReporter()
        mykrobereporter.add_input_files({'CSV': mykrobe.tool_outputs['CSV']})
        mykrobereporter.add_input_informs({'mykrobe': mykrobe.informs})
        mykrobereporter.run(self.running_dir)
        self.assertGreater(
            len(mykrobereporter.tool_outputs['HTML'][0].value.to_html()), 0
        )

        # Save the report
        html_out = self.running_dir / 'report.html'
        with html_out.open('w') as handle:
            handle.write(mykrobereporter.tool_outputs['HTML'][0].value.to_html())
        logging.info(f'Output report created: {html_out}')


if __name__ == '__main__':
    unittest.main()
