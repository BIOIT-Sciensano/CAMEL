import json
import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.scripts.sccmecfinder.mainsccmecfinder import MainSCCmecFinder


class TestSCCmecFinder(CamelTestSuite):
    """
    Tests the SCCmecFinder main script.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('sccmec')
    input_fasta_file = test_file_dir / 'MSSA476.fasta'

    def test_sccmecfinder_fasta(self) -> None:
        """
        Tests the SCCmecFinder tool with a FASTA input file.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(TestSCCmecFinder.input_fasta_file),
            '--input-type', 'fasta',
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--db-mec-genes', '/db/gene_detection/SCCmec_genes',
            '--profiles-mec-genes', '/db/pipelines/saureus/SCCmec_genes/profiles.yml'
        ]
        main = MainSCCmecFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_sccmecfinder_fasta_json_out(self) -> None:
        """
        Tests the SCCmecFinder tool with an input assembly and JSON output.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        output_json = self.running_dir / 'report' / 'informs.json'
        args = [
            '--fasta', str(TestSCCmecFinder.input_fasta_file),
            '--input-type', 'fasta',
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--output-json', str(output_json),
            '--working-dir', str(self.running_dir),
            '--db-mec-genes', '/db/gene_detection/SCCmec_genes',
            '--profiles-mec-genes', '/db/pipelines/saureus/SCCmec_genes/profiles.yml'
        ]
        main = MainSCCmecFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)
        with output_json.open() as handle:
            informs = json.load(handle)
            self.assertGreater(len(informs), 0)


if __name__ == '__main__':
    unittest.main()
