import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.tools.abritamr.abritamrreport import AbriTAMRReport
from camel.app.tools.abritamr.abritamrreporter import AbriTAMRReporter
from camel.app.tools.abritamr.abritamrrun import AbriTAMRRun
from camel.scripts.abritamr.mainabritamr import MainAbriTAMR


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
        abritamr = AbriTAMRRun(self.camel)
        abritamr.add_input_files({
            'FASTA': [ToolIOFile(Path(self.input_fasta_file))],
            'DIR_AMRF': [ToolIODirectory(Path('/db/abritamr/amrfinderplus_data'))]
        })
        abritamr.update_parameters(species='Salmonella')
        abritamr.run(self.running_dir)
        self.verify_output_files(abritamr, 'TXT_MATCHES')
        self.verify_output_files(abritamr, 'TXT_PARTIALS')

    def test_abritamr_report_and_reporter(self) -> None:
        """
        Tests basic AbriTAMR report and reporter.
        (abritamr report is a commandline tool and therefore not the camel reporter.)
        :return: None
        """
        # test abritamr report
        abritamr = AbriTAMRReport(self.camel)
        abritamr.add_input_files({
            'TXT_MDU_QC': [ToolIOFile(self.qc_file)],
            'TXT_MATCHES': [ToolIOFile(self.summary_matches)],
            'TXT_PARTIALS': [ToolIOFile(self.summary_partials)]
        })
        informs_abritamr = {'_name': 'test', 'species': 'Salmonella'}
        SnakemakeUtils.dump_object(informs_abritamr, Path('./informs.dummy'))
        abritamr.add_input_informs({'ABRITAMR_RUN': informs_abritamr})
        abritamr.run(self.running_dir)
        self.verify_output_files(abritamr, 'REPORT_ABRITAMR')

        # test abritamr reporter
        abritamrreporter = AbriTAMRReporter(self.camel)
        abritamrreporter.add_input_files({
            'TXT_MATCHES': [ToolIOFile(self.summary_matches)],
            'TXT_PARTIALS': [ToolIOFile(self.summary_partials)],
            'REPORT_ABRITAMR':  abritamr.tool_outputs['REPORT_ABRITAMR']
        })
        abritamrreporter.add_input_informs({'ABRITAMR_RUN': informs_abritamr})
        abritamrreporter.run(self.running_dir)
        output_section = abritamrreporter.tool_outputs['VAL_HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)
        CamelTestSuite.export_report_section(output_section, self.running_dir / 'report')

    def test_abritamr_standalone(self) -> None:
        """
        Tests the AbriTAMR standalone pipeline with fasta files.
        :return: None
        """
        path_report_html = self.running_dir / 'out' / 'report.html'
        path_report_tsv = self.running_dir / 'out' / 'report.tsv'

        args = [
            '--fasta', str(self.input_fasta_file),
            '--output-html', str(path_report_html),
            '--output-dir', str(path_report_html.parent),
            '--working-dir', str(self.running_dir),
            '--output-tsv', str(path_report_tsv),
            '--input-type', 'fasta',
            '--threads', '2',
            '--species', 'Salmonella'
        ]
        main = MainAbriTAMR(args)
        main.run()
        self.assertGreater(path_report_html.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
