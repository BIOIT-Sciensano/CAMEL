import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.tools.spotyping.spotyping import SpoTyping
from camel.app.tools.spotyping.spotypingreporter import SpoTypingReporter


class TestSpoTyping(CamelTestSuite):
    """
    Tests the SpoTyping tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir()
    input_fastq_pe = [
        test_file_dir / 'pipelines' / 'Myco-DRR041783-ds_1.fastq.gz',
        test_file_dir / 'pipelines' / 'Myco-DRR041783-ds_2.fastq.gz'
    ]
    input_fasta = test_file_dir / 'pipelines' / 'Myco-DRR041783-ds.fasta'

    def test_spotyping(self) -> None:
        """
        Tests the SpoTyping tool.
        :return: None
        """
        spotyping = SpoTyping()
        spotyping.add_input_files({
            'FASTQ': [ToolIOFile(file_) for file_ in self.input_fastq_pe],
            'INFORMS_spoligo_param': [ToolIOValue({
                'min_strict': 3,
                'min_relaxed': 3,
                'downsample_factor': 'NA'
            })]
        })
        spotyping.run(self.running_dir)
        self.assertIn('VAL_type_binary', spotyping.tool_outputs)
        self.verify_output_files(spotyping, 'LOG')

    def test_spotyping_swift_off(self) -> None:
        """
        Tests the SpoTyping tool with the 'swift' parameter disabled.
        :return: None
        """
        spotyping = SpoTyping()
        spotyping.add_input_files({
            'FASTQ': [ToolIOFile(file_) for file_ in self.input_fastq_pe],
            'INFORMS_spoligo_param': [ToolIOValue({
                'min_strict': 3,
                'min_relaxed': 3,
                'downsample_factor': 'NA'
            })]
        })
        spotyping.update_parameters(swift='off')
        spotyping.run(self.running_dir)
        self.assertIn('VAL_type_binary', spotyping.tool_outputs)
        self.verify_output_files(spotyping, 'LOG')

    def test_spotyping_fasta(self) -> None:
        """
        Tests the SpoTyping tool with FASTA input.
        :return: None
        """
        spotyping = SpoTyping()
        spotyping.add_input_files({
            'FASTA': [ToolIOFile(self.input_fasta)]
        })
        spotyping.update_parameters(fasta=None)
        spotyping.run(self.running_dir)
        self.assertIn('VAL_type_binary', spotyping.tool_outputs)
        self.verify_output_files(spotyping, 'LOG')

    def test_spotyping_reporter(self) -> None:
        """
        Tests the SpoTyping tool and reporter
        :return None
        """
        # Run the tool
        spotyping = SpoTyping()
        informs_spoligo_param = [ToolIOValue({
                'min_strict': 3,
                'min_relaxed': 3,
                'downsample_factor': 'NA'
            })]
        spotyping.add_input_files({
            'FASTQ': [ToolIOFile(file_) for file_ in self.input_fastq_pe],
            'INFORMS_spoligo_param': informs_spoligo_param
        })
        spotyping.run(self.running_dir)

        # Run the reporter
        reporter = SpoTypingReporter()
        reporter.add_input_files({
            'VAL_type_binary': spotyping.tool_outputs['VAL_type_binary'],
            'VAL_type_octal': spotyping.tool_outputs['VAL_type_octal'],
        })
        reporter.add_input_informs({
            'spotyping': spotyping.informs,
            'spoligo_param': informs_spoligo_param
        })
        reporter.run(self.running_dir)
        self.assertGreater(len(reporter.tool_outputs['VAL_HTML'][0].value.to_html()), 0)

    def test_spotyping_reporter_fasta(self) -> None:
        """
        Tests the SpoTyping tool and reporter with FASTA input
        :return None
        """
        # Run the tool
        spotyping = SpoTyping()
        spotyping.add_input_files({
            'FASTA': [ToolIOFile(self.input_fasta)]
        })
        spotyping.update_parameters(fasta=None)
        spotyping.run(self.running_dir)

        # Run the reporter
        reporter = SpoTypingReporter()
        reporter.add_input_files({
            'VAL_type_binary': spotyping.tool_outputs['VAL_type_binary'],
            'VAL_type_octal': spotyping.tool_outputs['VAL_type_octal'],
        })
        reporter.add_input_informs({'spotyping': spotyping.informs})
        reporter.run(self.running_dir)
        self.assertGreater(len(reporter.tool_outputs['VAL_HTML'][0].value.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
