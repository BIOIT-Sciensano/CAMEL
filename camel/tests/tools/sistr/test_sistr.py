import unittest
from pathlib import Path

from camel.app.config import config
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.pipelines.salmonella.sistr import Sistr
from camel.app.tools.pipelines.salmonella.sistrreporter import SISTRReporter


class TestSistr(CamelTestSuite):
    """
    Tests the Sistr tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('salmonella')
    input_fasta_file = test_file_dir / 'assembly_filtered.fasta'
    db_path = Path(config.dir_db / 'SISTR/1.1.1/data')

    def test_sistr(self) -> None:
        """
        Tests basic Sistr run.
        :return: None
        """
        sistr_tool = Sistr()
        sistr_tool.add_input_files({
            'FASTA': [ToolIOFile(Path(TestSistr.input_fasta_file))],
            'DIR': [ToolIODirectory(self.db_path)]
        })
        sistr_tool.run(self.running_dir)
        self.verify_output_files(sistr_tool, 'JSON')
        self.assertIn('_name', sistr_tool.informs)

    def test_sistr_reporter(self) -> None:
        """
        Tests Sistr reporter.
        :return: None
        """
        sistr_tool = Sistr()
        sistr_tool.add_input_files({
            'FASTA': [ToolIOFile(Path(TestSistr.input_fasta_file))],
            'DIR': [ToolIODirectory(self.db_path)]
        })
        sistr_tool.run(self.running_dir)
        self.verify_output_files(sistr_tool, 'JSON')

        sistr_reporter = SISTRReporter()
        sistr_reporter.add_input_files({
            'JSON_SISTR': sistr_tool.tool_outputs['JSON'],
            'DIR_sistr': [ToolIODirectory(self.db_path)]})
        sistr_reporter.add_input_informs({'serotyping_sistr': sistr_tool.informs})
        sistr_reporter.run(self.running_dir)
        output_section = sistr_reporter.tool_outputs['VAL_HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
