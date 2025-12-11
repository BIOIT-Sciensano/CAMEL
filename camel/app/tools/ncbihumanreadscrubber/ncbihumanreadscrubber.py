import re
import tempfile
from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.config import config
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.core.tool import Tool


class NcbiHumanReadScrubber(Tool):
    """
    NCBI human read scrubbing tool, also called HRRT or human read removal tool.
    """

    def __init__(self) -> None:
        """
        Initializes the HRRT.
        """
        super().__init__('HRRT', '2.2.1')

    def _execute_tool(self) -> None:
        """
        Runs the HRRT.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='hrrt_', dir=config.dir_temp) as dir_temp:
            self.__build_command(Path(dir_temp))
            self._execute_command(env={'TMPDIR': dir_temp})
            self._parse_stderr()
            self.__set_output()

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        if 'FASTQ_SE' not in self._tool_inputs or len(self._tool_inputs['FASTQ_SE']) != 1:
            raise InvalidToolInputError("Required FASTQ_SE input file is missing for human read scrubber.")
        super()._check_input()

    def __build_command(self, dir_temp: Path) -> None:
        """
        Builds the command line call to execute HRRT.
        Export_human_reads and outputfile_removed linked, adds args -r -u and path if export_human_reads = true
        :param dir_temp: path to the temporary directory
        :return: None
        """
        parts = [
            self._tool_command,
            *self._build_options(excluded_parameters=['interleaved', 'export_human_reads', 'outputfile_removed']),
            self._parameters['interleaved'].option if 'interleaved' in self._parameters else '',
            '-i', str(self._tool_inputs['FASTQ_SE'][0].path)
        ]
        if 'DB' in self._tool_inputs:
            parts.extend(['-d', str(self._tool_inputs['DB'][0].path)])
        if 'export_human_reads' in self._parameters:
            parts.extend([
                self._parameters['export_human_reads'].option,
                self._parameters['outputfile_removed'].option,
                str(Path(self._folder, self._parameters['outputfile_removed'].value))
            ])
        self._command.command = ' '.join(parts)

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def __set_output(self) -> None:
        """
        Set the output of HRRT.
        :return: None
        """
        # Scrubbed reads
        path_out = self.folder / self._parameters['outputfile'].value
        self._tool_outputs['FASTQ_SCRUBBED'] = [ToolIOFile(path_out)]

        # Human reads
        if 'export_human_reads' in self._parameters:
            if self._informs.get('statistics').get('count_removed') == 0:
                logger.warning('Human read export enabled, but no human reads found')
                self._tool_outputs['FASTQ_REMOVED'] = []
            else:
                path_removed = self.folder / self._parameters['outputfile_removed'].value
                self._tool_outputs['FASTQ_REMOVED'] = [ToolIOFile(path_removed)]

    def _parse_stderr(self) -> None:
        """
        Parses the command's stderr to retrieve the statistics about how many reads/contigs were removed.
        :return: None
        """
        count_removed = None
        count_total = None
        for line in self._command.stderr.splitlines():
            # Define the regular expression pattern
            pattern_reads_removed = r'^(\d+)\s+spot\(s\) masked or removed\.$'
            pattern_reads_total = r'total spot count: (\d+)'

            # Try to match the pattern in the current line
            if count_removed is None and re.match(pattern_reads_removed, line):
                # Extract the matched integer
                count_removed = int((re.match(pattern_reads_removed, line)).group(1))
            elif count_total is None and re.search(pattern_reads_total, line):
                # Extract the matched integer
                count_total = int((re.search(pattern_reads_total, line)).group(1))
            elif count_removed is not None and count_total is not None:
                break
        if count_removed is None or count_total is None:
            raise ValueError("The statistics of the human read scrubbing could not be obtained.")
        else:
            self._informs['statistics'] = {'count_removed': count_removed, 'count_total': count_total}
