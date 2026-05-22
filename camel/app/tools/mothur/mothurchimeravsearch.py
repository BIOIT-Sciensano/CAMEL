import re

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.loggers import logger
from camel.app.tools.mothur.mothur import Mothur

# When chimera.vsearch is run with a grouped count table, mothur tries to open per-sample
# split files (e.g. file.A.fasta / file.A.count_table) that don't exist separately.
# These [ERROR] lines are non-fatal — mothur continues and produces correct output —
# but they cause mothur to exit with a non-zero code, which would otherwise be treated
# as a tool failure.
_SPLIT_FILE_ERROR = re.compile(r'\[ERROR]: Could not open .+\.\w+\.(fasta|count_table)')


class MothurChimeraVsearch(Mothur):
    """
    The chimera.vsearch command reads a fasta file and reference file or a fasta and name or count file and outputs
    potentially chimeric sequences. The vsearch program is donated to the public domain,
    https://github.com/torognes/vsearch
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_chimera_vsearch')
        self._required_input = ['FASTA']
        self._optional_input = ['FASTA_Ref', 'TSV_Groups']

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid cfr. implemented in the superclass.
        Additionally:
        - Either TSV_Names or TSV_Counts is required
        - The use of TSV_Names is not yet implemented (lack of documentation)
        :return: None
        """
        if 'TSV_Counts' in self._tool_inputs:
            self._required_input.append('TSV_Counts')
        elif 'TSV_Names' in self._tool_inputs:
            self._required_input.append('TSV_Names')
        else:
            raise InvalidToolInputError('Either TSV_Counts or TSV_Names is required')
        super()._check_input()

    def _check_command_output(self, command: Command) -> None:
        """
        Checks command output, tolerating the non-fatal per-sample split-file errors that
        chimera.vsearch emits when a grouped count table is used but the per-sample FASTA /
        count_table files do not exist separately on disk. All other non-zero exits raise
        a ToolExecutionError.
        :return: None
        """
        if command.exit_code == 0:
            return
        error_lines = [line for line in command.stdout.splitlines() if line.startswith('[ERROR]')]
        if error_lines and all(_SPLIT_FILE_ERROR.search(line) for line in error_lines):
            for line in error_lines:
                logger.warning(f'Ignoring non-fatal per-sample split-file error: {line.strip()}')
            return
        raise ToolExecutionError(self.name, f"Error executing '{self.name}', exit code: {command.exit_code}")

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = [f"fasta={self._tool_inputs['FASTA'][0]}"]
        if 'TSV_Counts' in self._tool_inputs:
            items.append(f"count={self._tool_inputs['TSV_Counts'][0]}")
        elif 'TSV_Names' in self._tool_inputs:
            items.append(f"name={self._tool_inputs['TSV_Names'][0]}")
        if 'TSV_Groups' in self._tool_inputs:
            items.append(f"group={self._tool_inputs['TSV_Groups'][0]}")
        if 'FASTA_Ref' in self._tool_inputs:
            items.append(f"reference={self._tool_inputs['FASTA_Ref'][0]}")
        items.append(f"outputdir={self._folder}")
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them.
        chimera.vsearch always writes the non-chimeric sequences to .denovo.vsearch.fasta and
        .denovo.vsearch.count_table directly in outputdir, regardless of the dereplicate setting.
        When dereplicate=false mothur also runs remove.seqs internally which produces .pick.*
        intermediates, but those are not guaranteed to land in outputdir and must not be used.
        :return: None
        """
        basename = self._get_basename()
        self._tool_outputs['TSV_Chimeras'] = [ToolIOFile(basename.with_suffix('.denovo.vsearch.chimeras'))]
        self._tool_outputs['TSV_Accnos'] = [ToolIOFile(basename.with_suffix('.denovo.vsearch.accnos'))]
        self._tool_outputs['FASTA'] = [ToolIOFile(basename.with_suffix('.denovo.vsearch.fasta'))]
        if 'TSV_Counts' in self._tool_inputs:
            self._tool_outputs['TSV_Counts'] = [ToolIOFile(basename.with_suffix('.denovo.vsearch.count_table'))]
        if 'TSV_Names' in self._tool_inputs:
            raise RuntimeError('The use of a names file is not yet implemented for chimera.vsearch as the outputs are unknown!')
