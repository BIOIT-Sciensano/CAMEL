import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.config import config
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.nextclade.nextclade import Nextclade
from camel.app.tools.nextclade.nextcladereporter import NextcladeReporter


class TestNextclade(CamelTestSuite):
    """
    Tests the Nextclade tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('nextclade')

    def test_nextclade_h1n1(self) -> None:
        """
        Tests the Nextclade tool on influenza H1N1 data.
        :return: None
        """
        # Run Nextclade
        nextclade = Nextclade()
        nextclade.add_input_files({
            'FASTA': [ToolIOFile(TestNextclade.test_file_dir / 'sequences_h1n1pdm.fasta')],
            'DB': [ToolIODirectory(config.dir_db / 'nextclade' / 'flu_h1n1pdm_ha')]
        })
        nextclade.run(self.running_dir)
        self.verify_output_files(nextclade, 'CSV')
        self.assertGreater(len(nextclade.informs['results']), 0)

        # Run reporter
        reporter = NextcladeReporter()
        reporter.add_input_files({
            'CSV': nextclade._tool_outputs['CSV'],
            'DB': [ToolIODirectory(config.dir_db / 'nextclade' / 'flu_h1n1pdm_ha')]
        })
        reporter.add_input_informs({'nextclade': nextclade.informs})
        reporter.run()
        self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)

    def test_nextclade_h3n2(self) -> None:
        """
        Tests the Nextclade tool on influenza H3N2 data.
        :return: None
        """
        nextclade = Nextclade()
        nextclade.add_input_files({
            'FASTA': [ToolIOFile(TestNextclade.test_file_dir / 'sequences_h3n2.fasta')],
            'DB': [ToolIODirectory(Path(config.dir_db, 'nextclade', 'flu_h3n2_ha'))]
        })
        nextclade.run(self.running_dir)
        self.verify_output_files(nextclade, 'CSV')
        self.assertGreater(len(nextclade.informs['results']), 0)

        # Run reporter
        reporter = NextcladeReporter()
        reporter.add_input_files({
            'CSV': nextclade._tool_outputs['CSV'],
            'DB': [ToolIODirectory(Path(config.dir_db, 'nextclade', 'flu_h3n2_ha'))]
        })
        reporter.add_input_informs({'nextclade': nextclade.informs})
        reporter.run()
        self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)

    def test_nextclade_sars_cov_2(self) -> None:
        """
        Tests the Nextclade tool on SARS-CoV-2 data.
        :return: None
        """
        nextclade = Nextclade()
        nextclade.add_input_files({
            'FASTA': [ToolIOFile(TestNextclade.test_file_dir / 'sequences_sars-cov-2.fasta')],
            'DB': [ToolIODirectory(Path(config.dir_db, 'nextclade', 'sars-cov-2'))]
        })
        nextclade.run(self.running_dir)
        self.verify_output_files(nextclade, 'CSV')
        self.assertGreater(len(nextclade.informs['results']), 0)

        # Run reporter
        reporter = NextcladeReporter()
        reporter.add_input_files({
            'CSV': nextclade._tool_outputs['CSV'],
            'DB': [ToolIODirectory(Path(config.dir_db, 'nextclade', 'sars-cov-2'))]
        })
        reporter.add_input_informs({'nextclade': nextclade.informs})
        reporter.run()
        self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
