from pathlib import Path

import pandas as pd

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class ConFindr(Tool):
    """
    This program is designed to find bacterial intra-species contamination in raw Illumina data. It does this by looking
    for multiple alleles of core, single copy genes.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
                :return: None
        """
        super().__init__('ConFindr', '0.8.2')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if not any(x in self._tool_inputs for x in ('FASTQ_SE', 'FASTQ_PE')):
            raise InvalidToolInputError('FASTQ_(SE/PE) input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        dir_input = self.__symlink_input_files()
        dir_out = Path(self._folder) / 'out'
        self._command.command = ' '.join([
            self._tool_command,
            '--input_directory', str(dir_input),
            '--output_name', str(dir_out)
        ] + self._build_options())
        self._execute_command()
        self.__set_output(dir_out)
        self._check_output_content()

    def __symlink_input_files(self) -> Path:
        """
        Symlinks the input files.
        :return: Path to input directory
        """
        dir_working = Path(self._folder)
        if 'FASTQ_SE' in self._tool_inputs:
            is_gzipped = FileSystemHelper.is_gzipped(self._tool_inputs['FASTQ_SE'][0].path)
            extension = 'fastq' if not is_gzipped else 'fastq.gz'
            (dir_working / f'sample_R1.{extension}').symlink_to(self._tool_inputs['FASTQ_SE'][0].path)
        elif 'FASTQ_PE' in self._tool_inputs:
            is_gzipped = FileSystemHelper.is_gzipped(self._tool_inputs['FASTQ_PE'][0].path)
            extension = 'fastq' if not is_gzipped else 'fastq.gz'
            (dir_working / f'sample_R1.{extension}').symlink_to(self._tool_inputs['FASTQ_PE'][0].path)
            (dir_working / f'sample_R2.{extension}').symlink_to(self._tool_inputs['FASTQ_PE'][1].path)
        return dir_working

    def __set_output(self, dir_out: Path) -> None:
        """
        Sets the tool output for this tool.
        :param dir_out: Output directory
        :return: None
        """
        csv_file = dir_out / 'confindr_report.csv'
        self._tool_outputs['CSV'] = [ToolIOFile(csv_file)]
        self._informs.update(pd.read_csv(csv_file).to_dict('records')[0])

    def _check_output_content(self) -> None:
        """
        Checks the tool output content for errors.
        :return: None
        """
        if self._informs['Genus'] == 'Error processing sample':
            self._informs['has_error'] = True  # set, but not used for now, perhaps in future
            # important, setting NumContamSNVs to None will cause the tool to be noted as skipped in the overview
            self._informs['NumContamSNVs'] = None

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
