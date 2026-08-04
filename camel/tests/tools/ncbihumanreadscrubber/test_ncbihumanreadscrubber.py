import unittest

from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import fastqutils

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.ncbihumanreadscrubber.ncbihumanreadscrubber import (
    NcbiHumanReadScrubber,
)
from camel.app.tools.ncbihumanreadscrubber.ncbihumanreadscrubberreporter import (
    NcbiHumanReadScrubberReporter,
)


class TestNcbiHumanReadScrubber(CamelTestSuite):
    """
    Tests the HRRT tool.
    _nh: files without human reads
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('ncbi_human_read_scrubber')

    path_fq_se = test_file_dir / 'Myco-DRR041783-ds_1_subset.fastq'
    path_fq_ont = test_file_dir / 'minion_reads.fastq'
    path_fq_ont_nh = test_file_dir / 'minion_no_human.fastq'
    files_list = [path_fq_se, path_fq_ont, path_fq_ont_nh]
    output = [True, True, False]

    def test_scrubber_fq_se(self) -> None:
        """
        Tests the scrubber on a not gzipped single end FASTQ file.
        :return: None
        """
        for file, out in zip(TestNcbiHumanReadScrubber.files_list, TestNcbiHumanReadScrubber.output):
            scrubber = NcbiHumanReadScrubber()
            scrubber.add_input_files({'FASTQ_SE': [ToolIOFile(file)]})
            scrubber.update_parameters(interleaved=False, export_human_reads=True, outputfile=self.running_dir / 'test_scrubber_output.fastq', outputfile_removed=self.running_dir / 'test_scrubber_reads_removed.fastq')
            # Updated in the ncbi.py as excluded
            scrubber.run(self.running_dir)
            self.verify_output_files(scrubber, 'FASTQ_SCRUBBED', 1)
            self.verify_output_files(scrubber, 'FASTQ_REMOVED', 1) if out else self.verify_output_files(scrubber, 'FASTQ_REMOVED', 0)

            # Check that the input file is larger than the output file
            self.assertGreater(
                fastqutils.count_reads(file),
                fastqutils.count_reads(scrubber.tool_outputs['FASTQ_SCRUBBED'][0].path)
            ) if out else self.assertEqual(
                fastqutils.count_reads(file),
                fastqutils.count_reads(scrubber.tool_outputs['FASTQ_SCRUBBED'][0].path)
            )

            # Check if the informs were added
            self.assertIn('statistics', scrubber.informs)

    def test_scrubber_reporter_fq_se(self) -> None:
        """
        Tests the NCBI human read scrubber reporter class on a not gzipped single end FASTQ file.
        :return: None
        """
        for file, out in zip(TestNcbiHumanReadScrubber.files_list, TestNcbiHumanReadScrubber.output):
            scrubber = NcbiHumanReadScrubber()
            scrubber.add_input_files({'FASTQ_SE': [ToolIOFile(file)]})
            scrubber.update_parameters(interleaved=False, export_human_reads=True, outputfile=self.running_dir / 'test_scrubber_output.fastq', outputfile_removed=self.running_dir / 'test_scrubber_reads_removed.fastq')
            # Updated in the ncbi.py as excluded
            scrubber.run(self.running_dir)

            reporter = NcbiHumanReadScrubberReporter()
            reporter.add_input_informs({'SCRUBBER': scrubber.informs})
            reporter.add_input_files({'REMOVED': scrubber.tool_outputs['FASTQ_REMOVED']})
            reporter.update_parameters(input_format='fastq_se')
            reporter.run(self.running_dir)
            output_section = reporter.tool_outputs['HTML'][0].value
            self.assertGreater(len(output_section.to_html()), 0)

    def test_scrubbing_default(self) -> None:
        """
        Tests the NCBI human read scrubbing standalone pipeline with the default args.
        :return: None
        """
        for fastq, val in zip(TestNcbiHumanReadScrubber.files_list, TestNcbiHumanReadScrubber.output):
            scrubber_def = NcbiHumanReadScrubber()
            scrubber_def.add_input_files({'FASTQ_SE': [ToolIOFile(fastq)]})
            scrubber_def.update_parameters(
                interleaved=False,
                export_human_reads=False,
                outputfile=self.running_dir / 'test_scrubber_output.fastq'
            )
            # Updated in the ncbi.py as excluded
            scrubber_def.run(self.running_dir)
            self.verify_output_files(scrubber_def, 'FASTQ_SCRUBBED', 1)

            # Check that the input file is larger than the output file
            self.assertEqual(
                fastqutils.count_reads(fastq),
                fastqutils.count_reads(scrubber_def.tool_outputs['FASTQ_SCRUBBED'][0].path)
            ) if not val else self.assertGreater(
                fastqutils.count_reads(fastq),
                fastqutils.count_reads(scrubber_def.tool_outputs['FASTQ_SCRUBBED'][0].path))

            # Check if the informs were added
            self.assertIn('statistics', scrubber_def.informs)


if __name__ == '__main__':
    unittest.main()
