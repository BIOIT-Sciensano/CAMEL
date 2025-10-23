import unittest
from pathlib import Path

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.config import config
from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.nextclade3.nextclade3 import Nextclade3
from camel.app.tools.nextclade3.nextclade3reporter import Nextclade3Reporter


class TestNextclade3(CamelTestSuite):
    """
    Tests the Nextclade3 tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('nextclade')
    dir_db = config.dir_db / 'nextclade3'

    def test_nextclade_sars_cov_2(self) -> None:
        """
        Tests the Nextclade tool on SARS-CoV-2 data.
        :return: None
        """
        nextclade = Nextclade3()
        nextclade.add_input_files({
            'FASTA': [ToolIOFile(TestNextclade3.test_file_dir / 'sequences_sars-cov-2.fasta')],
            'DB': [ToolIODirectory(Path(config.dir_db, 'nextclade3', 'sars-cov-2'))]
        })
        dir_ = self.running_dir / 'genome'
        dir_.mkdir(parents=True, exist_ok=True)
        nextclade.run(dir_)
        self.verify_output_files(nextclade, 'TSV')
        self.assertGreater(len(nextclade.informs['results']), 0)

        # Run reporter
        reporter = Nextclade3Reporter()
        reporter.add_input_files({
            'TSV': nextclade.tool_outputs['TSV'],
            'DB': [ToolIODirectory(Path(config.dir_db, 'nextclade3', 'sars-cov-2'))]
        })
        reporter.add_input_informs({'nextclade': [nextclade.informs]})
        reporter.update_parameters(name='test_sample')
        reporter.run()
        self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
