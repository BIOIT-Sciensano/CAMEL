import logging
import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.snakemake import snakemakeutils
from camel.app.tools.pipelines.neisseria.characterizeneisseriacapsule import CharacterizeNeisseriaCapsule
from camel.app.tools.pipelines.neisseria.characterizeneisseriacapsulereporter import \
    CharacterizeNeisseriaCapsuleReporter


class TestCharacterizeNeisseriaCapsule(CamelTestSuite):
    """
    Tests the CharacterizeNeisseriaCapsule tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('characterize_neisseria_capsule')
    fasta_in = test_file_dir / 'Neisseria_S16BD00092.fasta'

    def test_characterize_neisseria_capsule(self) -> None:
        """
        Tests the CharacterizeNeisseriaCapsule tool.
        :return: None
        """
        serogroup_tool = CharacterizeNeisseriaCapsule()
        serogroup_tool.add_input_files({'FASTA': [ToolIOFile(TestCharacterizeNeisseriaCapsule.fasta_in)]})
        serogroup_tool.update_parameters(threads=2)
        serogroup_tool.run(self.running_dir)
        logging.info(f'Successfully processed: {TestCharacterizeNeisseriaCapsule.fasta_in}')
        self.verify_output_files(serogroup_tool, 'TSV')

    def test_characterize_neisseria_capsule_reporter(self) -> None:
        """
        Tests the CharacterizeNeisseriaCapsuleReporter tool.
        :return: None
        """
        # Run the tool
        serogroup_tool = CharacterizeNeisseriaCapsule()
        serogroup_tool.add_input_files({'FASTA': [ToolIOFile(TestCharacterizeNeisseriaCapsule.fasta_in)]})
        serogroup_tool.update_parameters(threads=2)
        serogroup_tool.run(self.running_dir)

        # Run the reporter
        reporter = CharacterizeNeisseriaCapsuleReporter()
        reporter.add_input_files({
            'TSV': serogroup_tool.tool_outputs['TSV'],
            'JSON': serogroup_tool.tool_outputs['JSON']
        })
        reporter.add_input_informs({'detector': serogroup_tool.informs})
        reporter.run(self.running_dir)
        self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)

        # Save the report in a pickle
        path_html_io = self.running_dir / 'html.iob'
        snakemakeutils.dump_object(reporter.tool_outputs['HTML'], path_html_io)
        logging.info(f'Report pickle saved to: {path_html_io}')


if __name__ == '__main__':
    unittest.main()
