import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.fastani.fastani import FastANI
from camel.app.tools.fastani.fastanireporter import FastANIReporter


class TestFastANI(CamelTestSuite):
    """
    Tests the FastANI tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('fastani')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'escherichia.fasta')
    FILE_FASTA_QUERY = ToolIOFile(test_file_dir / 'shigella.fasta')
    FILE_FASTA_QUERY_LIST = ToolIOFile(test_file_dir / 'query_list.list')
    FILE_FASTA_REF_LIST = ToolIOFile(test_file_dir / 'reference_list.list')

    def test_fastani(self) -> None:
        """
        Tests the FastANI tool with default options.
        :return: None
        """
        fastani = FastANI()
        fastani.add_input_files({'FASTA_Q': [TestFastANI.FILE_FASTA_QUERY], 'FASTA_R': [TestFastANI.FILE_FASTA_REF]})
        fastani.run(self.running_dir)
        self.verify_output_files(fastani, 'TSV')

    def test_fastani_reporter(self) -> None:
        """
        Tests the FastANI reporter class.
        :return: None
        """
        # Run FastANI
        fastani = FastANI()
        fastani.add_input_files({'FASTA_Q': [TestFastANI.FILE_FASTA_QUERY], 'FASTA_R': [TestFastANI.FILE_FASTA_REF]})
        fastani.run(self.running_dir)
        self.verify_output_files(fastani, 'TSV')

        # Run the FastANI reporter
        reporter = FastANIReporter()
        reporter.add_input_files({'TSV': fastani.tool_outputs['TSV']})
        reporter.add_input_informs({'fastani': fastani.informs})
        reporter.run(self.running_dir)
        output_section = reporter.tool_outputs['HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)

    def test_fastani_reporter_lists(self) -> None:
        """
        Tests the FastANI reporter class with lists as input.
        :return: None
        """
        # Run FastANI
        fastani = FastANI()
        fastani.add_input_files({
            'TSV_FASTA_Q': [TestFastANI.FILE_FASTA_QUERY_LIST],
            'TSV_FASTA_R': [TestFastANI.FILE_FASTA_REF_LIST]})
        fastani.run(self.running_dir)
        self.verify_output_files(fastani, 'TSV')

        # Run the FastANI reporter
        reporter = FastANIReporter()
        reporter.add_input_files({'TSV': fastani.tool_outputs['TSV']})
        reporter.add_input_informs({'fastani': fastani.informs})
        reporter.run(self.running_dir)
        output_section = reporter.tool_outputs['HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
