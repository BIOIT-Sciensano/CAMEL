import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.pharokka.pharokka import Pharokka
from camel.app.tools.pharokka.pharokkamultiplotter import PharokkaMultiplotter
from camel.app.tools.pharokka.pharokkareporter import PharokkaReporter


class TestPharokka(CamelTestSuite):
    """
    Tests the Pharokka tool.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('pharokka')
    #fasta = test_file_dir / 'phages.fasta'
    fasta = test_file_dir / 'S_20_721_prophage.fasta'
    #fasta = test_file_dir / 'pBAD33.fasta'

    def test_pharokka(self) -> None:
        """
        Tests the Pharokka tool.
        :return: None
        """
        pharokka = Pharokka()
        pharokka.add_input_files({'FASTA': [ToolIOFile(Path(TestPharokka.fasta))]})
        pharokka.run(self.running_dir)
        self.verify_output_files(pharokka, 'GBK')
        self.verify_output_files(pharokka, 'TSV_STATS')
        self.verify_output_files(pharokka, 'TSV_CARD')
        self.verify_output_files(pharokka, 'TSV_VFDB')
        self.verify_output_files(pharokka, 'TSV_INPHARED')

    def test_pharokka_multiplotter(self) -> None:
        """
        Tests the Pharokka multiplotter.
        :return: None
        """
        # Run the tool
        pharokka = Pharokka()
        pharokka.add_input_files({'FASTA': [ToolIOFile(Path(TestPharokka.fasta))]})
        pharokka.run(self.running_dir)

        # Run the multiplotter
        multiplotter = PharokkaMultiplotter()
        multiplotter.add_input_files({'GBK': pharokka.tool_outputs['GBK']})
        multiplotter.run(self.running_dir)

    def test_pharokka_reporter(self) -> None:
        """
        Tests the PharokkaReporter tool.
        :return: None
        """
        # Run Pharokka
        pharokka = Pharokka()
        pharokka.add_input_files({'FASTA': [ToolIOFile(Path(TestPharokka.fasta))]})
        pharokka.run(self.running_dir)

        # Run the multiplotter
        multiplotter = PharokkaMultiplotter()
        multiplotter.add_input_files({'GBK': pharokka.tool_outputs['GBK']})
        multiplotter.run(self.running_dir)

        # Run the reporter
        reporter = PharokkaReporter()
        reporter.add_input_files({
            'GBK': pharokka.tool_outputs['GBK'],
            'TSV_STATS': pharokka.tool_outputs['TSV_STATS'],
            'TSV_CARD': pharokka.tool_outputs['TSV_CARD'],
            'TSV_VFDB': pharokka.tool_outputs['TSV_VFDB'],
            'TSV_INPHARED': pharokka.tool_outputs['TSV_INPHARED'],
            'PNG': multiplotter.tool_outputs['PNG']
        })
        reporter.add_input_informs({'pharokka': pharokka.informs})
        reporter.run(self.running_dir)
        self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)
        CamelTestSuite.export_report_section(reporter.tool_outputs['HTML'][0].value, self.running_dir / 'report.html')


if __name__ == '__main__':
    unittest.main()
