from pathlib import Path
from typing import Callable, Optional

from camel.app.core.utils import toolutils
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class FastpReporter(Tool):
    """
    Tool to create HTML reports for fastp.
    """

    COLUMNS_STATS = [
        {'key': 'q20_bases', 'name': 'Q20 bases', 'fmt': lambda x: f'{x:,}'},
        {'key': 'q20_rate', 'name': 'Q20 rate', 'fmt': lambda x: f'{x:.2f}'},
        {'key': 'q30_bases', 'name': 'Q30 bases', 'fmt': lambda x: f'{x:,}'},
        {'key': 'q30_rate', 'name': 'Q30 rate', 'fmt': lambda x: f'{x:.2f}'},
        {'key': 'read1_mean_length', 'name': 'Fwd. mean length', 'fmt': lambda x: f'{x:,}'},
        {'key': 'read2_mean_length', 'name': 'Rev. mean length', 'fmt': lambda x: f'{x:,}'},
        {'key': 'total_bases', 'name': 'Total bases', 'fmt': lambda x: f'{x:,}'},
        {'key': 'total_reads', 'name': 'Total reads', 'fmt': lambda x: f'{x:,}'}
    ]

    def __init__(self) -> None:
        """
        Initialize this tool.
        :return: None
        """
        super().__init__('fastp reporter', '0.1')
        self.__sub_folder = Path('read_trimming')
        self._report_section = None

    @staticmethod
    def _format_value(str_in: str, fmt: Optional[Callable] = None) -> str:
        """
        Formats the input value.
        :param str_in: input string
        :param fmt: formatting function
        :return: formatted string
        """
        if fmt is not None:
            return fmt(str_in)
        return str_in

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._report_section = HtmlReportSection('Read trimming', subtitle=self._input_informs['fastp']['_name'])
        self.__add_fastqc_reports('Pre-trimming', 'pre_trimming', 'HTML_pre')
        self.__add_fastp_section()
        self.__add_fastqc_reports('Post-trimming', 'post_trimming', 'HTML_post')
        if 'export_fastq' in self._parameters and self._parameters['export_fastq'].value is True:
            self.__add_trimmed_read_files()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section, False)]

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['HTML_pre', 'HTML_post', 'HTML'])
        if 'FASTQ_PE' not in self._tool_inputs and 'FASTQ_SE' not in self._tool_inputs:
            raise InvalidToolInputError("FASTQ input is required ('FASTQ_PE')")
        if 'fastp' not in self._input_informs:
            raise InvalidToolInputError("fastp inform are required ('fastp')")
        super()._check_input()

    def __add_fastqc_reports(self, header: str, dir_name: str, key: str) -> None:
        """
        Adds FastQC reports to the report. This function supports PE and SE input.
        :param header: Section header
        :param dir_name: Name of the directory to store files
        :param key: Input file key
        :return: None
        """
        self._report_section.add_header(header, 3)
        if len(self._tool_inputs[key]) == 1:
            relative_path = self.__sub_folder / dir_name / 'fastqc_report.html'
            self._report_section.add_file(self._tool_inputs[key][0].path, relative_path)
            self._report_section.add_link_to_file('FastQC report', relative_path)
        else:
            for fastqc_report, orientation in zip(self._tool_inputs[key], ('forward', 'reverse')):
                relative_path = self.__sub_folder / dir_name / f'fastqc_report_{orientation}.html'
                self._report_section.add_file(fastqc_report.path, relative_path)
                self._report_section.add_link_to_file(f'FastQC report ({orientation})', relative_path)

    def __add_fastp_section(self) -> None:
        """
        Adds the section with the fastp statistics.
        :return: None
        """
        self._report_section.add_header('Statistics', 3)
        header = ['Metric', 'Before filtering', 'After filtering']
        self._report_section.add_table([[
            key['name'],
            FastpReporter._format_value(self._input_informs['fastp']['summary']['before_filtering'][key['key']], key['fmt']),
            FastpReporter._format_value(self._input_informs['fastp']['summary']['after_filtering'][key['key']], key['fmt'])
        ] for key in FastpReporter.COLUMNS_STATS], header, [('class', 'data')])

        # Add a download link to the full report
        relative_path = self.__sub_folder / 'fastp.html'
        self._report_section.add_file(self._tool_inputs['HTML'][0].path, relative_path)
        self._report_section.add_link_to_file('fastp report (HTML)', relative_path)

    def __add_trimming_section_pe(self) -> None:
        """
        Adds the read trimming section for PE reads.
        :return: None
        """
        self._report_section.add_header('Read trimming', 3)
        table_data = [
            ['Input reads pairs:', FastpReporter.__reformat_inform(
                self._input_informs['trimming']['paired_reads_in'])],
            ['Both surviving:', FastpReporter.__reformat_inform(
                self._input_informs['trimming']['paired_reads_out'])],
            ['Forward only surviving:', FastpReporter.__reformat_inform(
                self._input_informs['trimming']['forward_only_reads'])],
            ['Reverse only surviving:', FastpReporter.__reformat_inform(
                self._input_informs['trimming']['reverse_only_reads'])],
            ['Dropped:', FastpReporter.__reformat_inform(
                self._input_informs['trimming']['reads_drop'])]
        ]
        self._report_section.add_table(table_data, table_attributes=[('class', 'information')])

    @staticmethod
    def __reformat_inform(input_str: str) -> str:
        """
        This function is used to reformat an inform value to a more readable format.
        This function also works when the percentage is omitted.
        E.g. 5241241 (10.02%) -> 5.241.241 (10.02%)
        :param input_str: Input string
        :return: Reformatted inform
        """
        parts = input_str.split(' ')
        if len(parts) == 1:
            return f'{int(parts[0]):,}'
        elif len(parts) == 2:
            return f'{int(parts[0]):,} {parts[1]}'
        raise ValueError(f"Cannot parse: {input_str}")

    def __add_trimmed_read_files(self) -> None:
        """
        Saves the PE trimmed reads in the output directory and adds them to the report.
        :return: None
        """
        self._report_section.add_header('Downloads', 3)
        for trimmed_reads_paired, orientation in zip(self._tool_inputs['FASTQ_PE'], ('forward', 'reverse')):
            relative_path = self.__sub_folder / 'read_trimming' / trimmed_reads_paired.path.name
            self._report_section.add_file(trimmed_reads_paired.path, relative_path)
            self._report_section.add_link_to_file(f'Trimmed {orientation} paired reads (FASTQ.GZ)', relative_path)
        self._report_section.add_line_break()

        if 'FASTQ_SE_FWD' in self._tool_inputs and len(self._tool_inputs['FASTQ_SE_FWD']) > 0:
            relative_path = self.__sub_folder / 'read_trimming' / self._tool_inputs['FASTQ_SE_FWD'][0].path.name
            self._report_section.add_file(self._tool_inputs['FASTQ_SE_FWD'][0].path, relative_path)
            self._report_section.add_link_to_file('Trimmed unpaired forwards reads (FASTQ.GZ)', relative_path)
        else:
            self._report_section.add_link_to_file('Trimmed unpaired forwards reads', None)

        if 'FASTQ_SE_REV' in self._tool_inputs and len(self._tool_inputs['FASTQ_SE_REV']) > 0:
            relative_path = self.__sub_folder / 'read_trimming' / self._tool_inputs['FASTQ_SE_REV'][0].path.name
            self._report_section.add_file(self._tool_inputs['FASTQ_SE_REV'][0].path, relative_path)
            self._report_section.add_link_to_file('Trimmed unpaired reverse reads (FASTQ.GZ)', relative_path)
        else:
            self._report_section.add_link_to_file('Trimmed unpaired reverse reads', None)
