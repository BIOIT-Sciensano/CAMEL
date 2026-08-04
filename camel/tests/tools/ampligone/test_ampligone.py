import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.config import config
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.ampligone.ampligone import AmpliGone
from camel.app.tools.ampligone.ampligonefasta2bed import AmpliGoneFasta2Bed
from camel.app.tools.ampligone.ampligonefasta2bedreporter import (
    AmpliGoneFasta2BedReporter,
)
from camel.app.tools.ampligone.ampligonereporter import AmpliGoneReporter


class TestAmpliGone(CamelTestSuite):
    """
    Tests the AmpliGone tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines', 'viral_consensus')
    FILE_FASTQ = ToolIOFile(test_file_dir / 'ESIB_EQA_2023.SARS1.01.fastq.gz')
    FILE_FASTA_PRIMERS = ToolIOFile(test_file_dir / 'ESIB_EQA_2023.SARS1.primers.fasta')
    FILE_FASTA_REF = ToolIOFile(Path(
        config.dir_db, 'pipelines', 'viral_consensus', 'ref_genomes', 'sars_cov_2-Wuhan-Hu-1.fasta'))

    def test_ampligone_ont(self) -> None:
        """
        Tests the AmpliGone tool with default options and ONT data.
        :return: None
        """
        ampligone = AmpliGone()
        ampligone.add_input_files({
            'FASTA_ref': [TestAmpliGone.FILE_FASTA_REF],
            'FASTA_primers': [TestAmpliGone.FILE_FASTA_PRIMERS],
            'FASTQ': [TestAmpliGone.FILE_FASTQ],
        })
        ampligone.update_parameters(amplicon_type='fragmented', error_rate='0.1')
        ampligone.run(self.running_dir)
        self.verify_output_files(ampligone, 'FASTQ')
        self.verify_output_files(ampligone, 'BED')
        self.assertIn('nucleotides_removed', ampligone.informs)
        self.assertIn('percentage_removed', ampligone.informs)

    def test_ampligone_reporter_ont(self) -> None:
        """
        Tests the AmpliGone tool with default options and ONT data.
        :return: None
        """
        # Run AmpliGone
        ampligone = AmpliGone()
        ampligone.add_input_files({
            'FASTA_ref': [TestAmpliGone.FILE_FASTA_REF],
            'FASTA_primers': [TestAmpliGone.FILE_FASTA_PRIMERS],
            'FASTQ': [TestAmpliGone.FILE_FASTQ],
        })
        ampligone.update_parameters(amplicon_type='fragmented', error_rate='0.1')
        ampligone.run(self.running_dir)

        # Run reporter
        reporter = AmpliGoneReporter()
        reporter.add_input_files({'BED': ampligone.tool_outputs['BED']})
        reporter.add_input_informs({'ampligone': ampligone.informs})
        reporter.run(self.running_dir)
        self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)

    def test_ampligone_fasta2bed(self) -> None:
        """
        Tests the AmpliGone BED file extraction of primer coordinates.
        :return: None
        """
        ampligone_fasta2bed = AmpliGoneFasta2Bed()
        ampligone_fasta2bed.add_input_files({
            'FASTA_ref': [TestAmpliGone.FILE_FASTA_REF],
            'FASTA_primers': [TestAmpliGone.FILE_FASTA_PRIMERS]
        })
        ampligone_fasta2bed.update_parameters(primer_mismatch_rate='0.1')
        ampligone_fasta2bed.run(self.running_dir)
        self.verify_output_files(ampligone_fasta2bed, 'BED')
        self.assertIn('primers_in', ampligone_fasta2bed.informs)
        self.assertIn('primers_out', ampligone_fasta2bed.informs)

    def test_ampligone_fasta2bed_reporter(self) -> None:
        """
        Tests the AmpliGone BED file extraction of primer coordinates.
        :return: None
        """
        ampligone_fasta2bed = AmpliGoneFasta2Bed()
        ampligone_fasta2bed.add_input_files({
            'FASTA_ref': [TestAmpliGone.FILE_FASTA_REF],
            'FASTA_primers': [TestAmpliGone.FILE_FASTA_PRIMERS]
        })
        ampligone_fasta2bed.update_parameters(primer_mismatch_rate='0.1')
        ampligone_fasta2bed.run(self.running_dir)

        # Run reporter
        reporter = AmpliGoneFasta2BedReporter()
        reporter.add_input_files({'BED': ampligone_fasta2bed.tool_outputs['BED']})
        reporter.add_input_informs({'ampligone': ampligone_fasta2bed.informs})
        reporter.run(self.running_dir)
        self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
