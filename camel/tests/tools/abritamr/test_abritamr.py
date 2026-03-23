import unittest
from pathlib import Path

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.abritamr.abritamrreport import AbriTAMRReport
from camel.app.tools.abritamr.abritamrreporter import AbriTAMRReporter
from camel.app.tools.abritamr.abritamrrun import AbriTAMRRun


class TestAbriTAMR(CamelTestSuite):
    """
    Tests the three AbriTAMR components; the tool called run, the tool called report, and the reporter.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('salmonella')
    qc_file = test_file_dir / 'qc_file.txt'
    summary_matches = test_file_dir / 'summary_matches.txt'
    summary_partials = test_file_dir / 'summary_partials.txt'
    input_fasta_file = test_file_dir / 'assembly_filtered.fasta'

    def test_abritamr_run(self) -> None:
        """
        Tests basic AbriTAMR run.
        :return: None
        """
        abritamr = AbriTAMRRun()
        abritamr.add_input_files({
            'FASTA': [ToolIOFile(Path(self.input_fasta_file))],
            'DIR_AMRF': [ToolIODirectory(Path('/db/amrfinder/v4/latest'))]
        })
        abritamr.update_parameters(species='Salmonella')
        abritamr.run(self.running_dir)
        self.verify_output_files(abritamr, 'TXT_matches')
        self.verify_output_files(abritamr, 'TXT_partials')

    def test_abritamr_report_and_reporter(self) -> None:
        """
        Tests basic AbriTAMR report and reporter.
        (abritamr report is a commandline tool and therefore not the camel reporter.)
        :return: None
        """
        # test abritamr report
        abritamr = AbriTAMRReport()
        abritamr.add_input_files({
            'TXT_mdu_qc': [ToolIOFile(self.qc_file)],
            'TXT_matches': [ToolIOFile(self.summary_matches)],
            'TXT_partials': [ToolIOFile(self.summary_partials)]
        })
        informs_abritamr = {'_name_full': 'test', 'species': 'Salmonella'}
        abritamr.add_input_informs({'abritamr_run': informs_abritamr})
        abritamr.run(self.running_dir)
        self.verify_output_files(abritamr, 'REPORT_abritamr')

        # test abritamr reporter
        abritamrreporter = AbriTAMRReporter()
        abritamrreporter.add_input_files({
            'TXT_matches': [ToolIOFile(self.summary_matches)],
            'TXT_partials': [ToolIOFile(self.summary_partials)],
            'REPORT_abritamr': abritamr.tool_outputs['REPORT_abritamr']
        })
        abritamrreporter.add_input_informs({'abritamr_run': informs_abritamr})
        abritamrreporter.run(self.running_dir)
        output_section = abritamrreporter.tool_outputs['VAL_HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)
        CamelTestSuite.export_report_section(output_section, self.running_dir / 'report')


if __name__ == '__main__':
    unittest.main()
