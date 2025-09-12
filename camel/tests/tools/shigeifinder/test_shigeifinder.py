import logging
import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.tools.pipelines.shigella.shigeifinder import ShigEiFinder
from camel.app.tools.pipelines.shigella.shigeifinderreporter import ShigEiFinderReporter


class TestShigEiFinder(CamelTestSuite):
    """
    Tests the ShigEiFinder tool.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('shigeifinder')
    fasta_in = test_file_dir / 'Shigella-S17BD07654.fasta'

    def test_shigeifinder(self) -> None:
        """
        Tests the ShigEiFinder tool.
        :return: None
        """
        serotype_tool = ShigEiFinder()
        serotype_tool.add_input_files({'FASTA': [ToolIOFile(TestShigEiFinder.fasta_in)]})
        serotype_tool.run(self.running_dir)
        logging.info(f'Successfully processed: {TestShigEiFinder.fasta_in}')
        self.verify_output_files(serotype_tool, 'TSV')

    def test_shigeifinder_reporter(self) -> None:
        """
        Tests the ShigEiFinderReporter tool.
        :return: None
        """
        # Run the tool
        serotype_tool = ShigEiFinder()
        serotype_tool.add_input_files({'FASTA': [ToolIOFile(TestShigEiFinder.fasta_in)]})
        serotype_tool.run(self.running_dir)

        # Run the reporter
        reporter = ShigEiFinderReporter()
        reporter.add_input_files({'TSV': serotype_tool.tool_outputs['TSV']})
        reporter.add_input_informs({'shigeifinder': serotype_tool.informs})
        reporter.run(self.running_dir)
        self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)

        # Save the report in a pickle
        path_html_io = self.running_dir / 'html.iob'
        SnakemakeUtils.dump_object(reporter.tool_outputs['HTML'], path_html_io)
        logging.info(f'Report pickle saved to: {path_html_io}')


if __name__ == '__main__':
    unittest.main()
