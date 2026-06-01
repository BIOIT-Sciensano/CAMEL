from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import fastqutils, sambamutils

from camel.app.core import toolutils
from camel.app.tools.samtools.samtoolsbase import SamtoolsBase


class SamtoolsFastq(SamtoolsBase):
    """
    Converts a BAM to a FASTQ.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('samtools fastq')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['BAM'])
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            str(self._tool_inputs['BAM'][0].path),
            *self._build_options()
        ])
        self._execute_command()
        self._set_output()
        self._verify_reads()

    def _set_output(self) -> None:
        """
        Collects the tool outputs.
        :return: None
        """
        # Split paired end output
        if all(x in self._params for x in ['read_1', 'read_2']):
            self._tool_outputs['FASTQ_PE'] = [
                ToolIOFile(self._folder / self.get_param_value('read_1')),
                ToolIOFile(self._folder / self.get_param_value('read_2'))
            ]

        # Combined output
        if 'output' in self._params:
            self._tool_outputs['FASTQ'] = [ToolIOFile(self._folder / self.get_param_value('output'))]

        # Non-PE reads
        if 'read_other' in self._params:
            self._tool_outputs['FASTQ_OTHER'] = [ToolIOFile(self._folder / self.get_param_value('read_other'))]

        # Singletons
        if 'singletons' in self._params:
            self._tool_outputs['FASTQ_SINGLETON'] = [ToolIOFile(self._folder / self.get_param_value('singletons'))]

        if len(self._tool_outputs) == 0:
            raise ValueError('Invalid output parameters')

    def _verify_reads(self) -> None:
        """
        Checks if all reads were converted.
        :return: None
        """
        nb_reads_in = sambamutils.get_record_count(self._tool_inputs['BAM'][0].path, primary_only=True)
        nb_reads_out = sum(sum(fastqutils.count_reads(x.path) for x in io) for k, io in self._tool_outputs.items())
        if nb_reads_in != nb_reads_out:
            raise ValueError(f'Number of input reads ({nb_reads_in}) does not match number of output reads ({nb_reads_out})')
