import logging
import unittest
from pathlib import Path

from camel.app.config import config
from camel.app.core.reports.htmlreport import HtmlReport
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.pipelines.neisseria.mendevarreporter import MenDeVARReporter
from camel.app.tools.pipelines.neisseria.mendevar import MenDeVAR
from camel.resources import CSS_STYLE


class TestMendevar(CamelTestSuite):
    """
    Tests the local MenDeVAR implementation.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mendevar')
    DB = Path(config.dir_db, 'pipelines', 'neisseria', 'mendevar_DB.txt')

    def test_mendevar(self) -> None:
        """
        Tests the MenDeVAR tool.
        :return: None
        """
        for tsv_path in TestMendevar.test_file_dir.glob('*.tsv'):
            # Create the output directory
            output_dir = self.running_dir / f'{tsv_path.stem}'
            output_dir.mkdir(parents=True)
            logging.info(f'Saving output in: {output_dir}')

            # Run MenDeVAR
            mendevar = MenDeVAR()
            mendevar.add_input_files({
                'TSV': [ToolIOFile(tsv_path)],
                'DB': [ToolIOFile(TestMendevar.DB)]
            })
            mendevar.update_parameters(output_directory=str(output_dir))
            mendevar.run(self.running_dir)
            self.verify_output_files(mendevar, 'TSV')
            logging.info(f'Successfully processed: {tsv_path}')

    def test_mendevar_reporter(self) -> None:
        """
        Tests the MenDeVAR reporter.
        :return: None
        """
        for tsv_path in TestMendevar.test_file_dir.glob('*.tsv'):
            # Create the output directory
            output_dir = self.running_dir / f'{tsv_path.stem}'
            output_dir.mkdir(parents=True)
            logging.info(f'Saving output in: {output_dir}')

            # Run MenDeVAR
            mendevar = MenDeVAR()
            mendevar.add_input_files({
                'TSV': [ToolIOFile(tsv_path)],
                'DB': [ToolIOFile(TestMendevar.DB)]
            })
            mendevar.update_parameters(output_directory=str(output_dir))
            mendevar.run(self.running_dir)
            self.verify_output_files(mendevar, 'TSV')
            logging.info(f'Successfully processed: {tsv_path}')

            # Run the reporter
            reporter = MenDeVARReporter()
            reporter.add_input_files({'TSV': mendevar.tool_outputs['TSV']})
            reporter.add_input_informs({'mendevar': mendevar.informs})
            reporter.run(self.running_dir)
            self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)

            # Save the report
            html_out = output_dir / 'report.html'
            report = HtmlReport(html_out)
            report.initialize('MenDeVar', CSS_STYLE)
            report.add_html_object(reporter.tool_outputs['HTML'][0].value)
            report.save()
            logging.info(f'Output report created: {html_out}')


if __name__ == '__main__':
    unittest.main()
