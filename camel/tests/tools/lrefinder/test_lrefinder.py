from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.lrefinder.lrefinder import LREFinder
from camel.app.tools.lrefinder.lrefinderreporter import LREFinderReporter


class TestLREFinder(CamelTestSuite):
    """
    Tests the LREFinder tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('lrefinder')
    input_fastq_pe = [test_file_dir / 'SRR12362697-ds_1.fastq.gz', test_file_dir / 'SRR12362697-ds_2.fastq.gz']
    input_fastq_se = test_file_dir / 'SRR17731943_ont-ds.fastq.gz'

    def test_lrefinder(self) -> None:
        """
        Tests the LREFinder tool.
        :return: None
        """
        lrefinder = LREFinder()
        lrefinder.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in TestLREFinder.input_fastq_pe]
        })
        lrefinder.run(self.running_dir)
        for key in ['species', 'genes', 'mutations']:
            self.assertIn(key, lrefinder.informs)

    def test_lrefinder_se(self) -> None:
        """
        Tests the LREFinder tool with FASTQ_SE input (like ont)
        :return: None
        """
        lrefinder = LREFinder()
        lrefinder.add_input_files({
            'FASTQ_SE': [ToolIOFile(TestLREFinder.input_fastq_se)]
        })
        lrefinder.run(self.running_dir)
        for key in ['species', 'genes', 'mutations']:
            self.assertIn(key, lrefinder.informs)

    def test_lrefinder_with_report(self) -> None:
        """
        Tests the LREFinder tool with the reporter.
        :return: None
        """
        lrefinder = LREFinder()
        lrefinder.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in TestLREFinder.input_fastq_pe]
        })
        lrefinder.run(self.running_dir)
        for key in ['species', 'genes', 'mutations']:
            self.assertIn(key, lrefinder.informs)

        reporter = LREFinderReporter()
        reporter.add_input_informs({'lrefinder': lrefinder.informs})
        reporter.run(self.running_dir)
        self.assertIn('HTML', reporter.tool_outputs)
