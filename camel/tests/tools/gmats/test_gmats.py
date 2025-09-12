import logging
import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.pipelines.neisseria.gmatsreporter import GMatsReporter
from camel.app.tools.pipelines.neisseria.gmats import GMats


class TestGMats(CamelTestSuite):
    """
    Tests the gMATS tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('gmats')
    DB = Path(CamelTestSuite.camel.config['db_root'], 'pipelines', 'neisseria', 'gMATS_DB.txt')

    def test_gmats(self) -> None:
        """
        Tests the gMATS tool.
        :return: None
        """
        for tsv_path in TestGMats.test_file_dir.glob('*.tsv'):
            # Create the output directory
            output_dir = self.running_dir / f'{tsv_path.stem}'
            output_dir.mkdir(parents=True)
            logging.info(f'Saving output in: {output_dir}')

            # Run gMATS
            gmats = GMats()
            gmats.add_input_files({
                'TSV': [ToolIOFile(tsv_path)],
                'DB': [ToolIOFile(TestGMats.DB)]
            })
            gmats.update_parameters(output_directory=str(output_dir))
            gmats.run(self.running_dir)
            self.verify_output_files(gmats, 'TSV')
            logging.info(f'Successfully processed: {tsv_path}')

    def test_gmats_reporter(self) -> None:
        """
        Tests the gMATS reporter.
        :return: None
        """
        for tsv_path in TestGMats.test_file_dir.glob('*.tsv'):
            # Create the output directory
            output_dir = self.running_dir / f'{tsv_path.stem}'
            output_dir.mkdir(parents=True)
            logging.info(f'Saving output in: {output_dir}')

            # Run gMATS
            gmats = GMats()
            gmats.add_input_files({
                'TSV': [ToolIOFile(tsv_path)],
                'DB': [ToolIOFile(TestGMats.DB)]
            })
            gmats.update_parameters(output_directory=str(output_dir))
            gmats.run(self.running_dir)
            self.verify_output_files(gmats, 'TSV')
            logging.info(f'Successfully processed: {tsv_path}')

            # Run the reporter
            reporter = GMatsReporter()
            reporter.add_input_files({'TSV': gmats.tool_outputs['TSV']})
            reporter.add_input_informs({'gmats': gmats.informs})
            reporter.run(self.running_dir)
            self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)

            # Save the report
            html_out = output_dir / 'report.html'
            with html_out.open('w') as handle:
                handle.write(reporter.tool_outputs['HTML'][0].value.to_html())
            logging.info(f'Output report created: {html_out}')


if __name__ == '__main__':
    unittest.main()
