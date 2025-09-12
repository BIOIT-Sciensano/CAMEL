import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.pipelines.shigella.shigatyper import ShigaTyper
from camel.app.tools.pipelines.shigella.shigatyperreporter import ShigaTyperReporter
from camel.tests import minOSVersion


class TestShigaTyper(CamelTestSuite):
    """
    Tests the ShigaTyper tool and the related reporter tool.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines')
    fastq_fwd = test_file_dir / 'Shigella-S17BD07654_1.fastq.gz'
    fastq_rev = test_file_dir / 'Shigella-S17BD07654_2.fastq.gz'
    fastq_ont = test_file_dir / 'Shigella-SRR29782656_ont-ds.fastq.gz'

    @minOSVersion('jammy')
    def test_shigatyper(self) -> None:
        """
        Tests the ShigaTyper tool.
        :return: None
        """
        shigatyper = ShigaTyper()
        shigatyper.add_input_files({
            'FASTQ_PE': [ToolIOFile(Path(TestShigaTyper.fastq_fwd)), ToolIOFile(Path(TestShigaTyper.fastq_rev))]})
        shigatyper.run(self.running_dir)
        self.verify_output_files(shigatyper, 'TSV')
        self.verify_output_files(shigatyper, 'TSV_HITS')

    @minOSVersion('jammy')
    def test_shigatyper_ont_input(self) -> None:
        """
        Tests the ShigaTyper tool with ONT input.
        :return: None
        """
        shigatyper = ShigaTyper()
        shigatyper.add_input_files({'FASTQ_SE': [ToolIOFile(TestShigaTyper.fastq_ont)]})
        shigatyper.update_parameters(ont=True)
        shigatyper.run(self.running_dir)
        self.verify_output_files(shigatyper, 'TSV')
        self.verify_output_files(shigatyper, 'TSV_HITS')

    @minOSVersion('jammy')
    def test_shigatyper_reporter(self) -> None:
        """
        Tests the ShigaTyperReporter tool.
        :return: None
        """
        # Run the tool
        shigatyper = ShigaTyper()
        shigatyper.add_input_files({
            'FASTQ_PE': [ToolIOFile(Path(TestShigaTyper.fastq_fwd)), ToolIOFile(Path(TestShigaTyper.fastq_rev))]})
        shigatyper.run(self.running_dir)

        # Run the reporter
        reporter = ShigaTyperReporter()
        reporter.add_input_files({
            'TSV': shigatyper.tool_outputs['TSV'],
            'TSV_HITS': shigatyper.tool_outputs['TSV_HITS']
        })
        reporter.add_input_informs({'shigatyper': shigatyper.informs})
        reporter.run(self.running_dir)

        # Check the output
        self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)
        CamelTestSuite.export_report_section(reporter.tool_outputs['HTML'][0].value, self.running_dir / 'report')


if __name__ == '__main__':
    unittest.main()
