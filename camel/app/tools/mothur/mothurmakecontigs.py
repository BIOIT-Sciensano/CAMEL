import re

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurMakeContigs(Mothur):
    """
    The make.contigs command reads a forward fastq file and a reverse fastq file and outputs new fasta and report files.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_make_contigs', '1.39.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - Only one input key is allowed
        - Only FASTQ_PE, FASTA_PE, TSV_File and TSV_Oligos are allowed as keys
        - Only 2 files allowed for PE and 1 for TSV_File and TSV_Oligos files
        :return: None
        """
        super()._check_input()
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidToolInputError('Too many input keys given for Mothur make.contigs: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.items():
            if key not in ['FASTQ_PE', 'FASTA_PE', 'TSV_File', 'TSV_Oligos']:
                raise InvalidToolInputError('Invalid input key given for Mothur make.contigs: {!r}'.format(self._tool_inputs))
            if key in ['FASTQ_PE', 'FASTA_PE']:
                if len(input_files) != 2:
                    raise InvalidToolInputError('Invalid number of files given for Mothur \
                                                         make.contigs: {!r}'.format(self._tool_inputs))
            else:
                if len(input_files) != 1:
                    raise InvalidToolInputError('Invalid number of files given for Mothur \
                                                         make.contigs: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and input/output directories
        Example: ffastq=FileR1.fastq, rfastq=FileR2.fastq, inputdir=/test/data/input/,
        outputdir=/test/data/outputdir
        :return: String with the input parameters
        """
        items = []
        for key, input_files in self._tool_inputs.items():
            if key == 'FASTQ_PE':
                items.append('ffastq={}, rfastq={}'.format(input_files[0], input_files[1]))
            elif key == 'FASTA_PE':
                items.append('ffasta={}, rfasta={}'.format(input_files[0], input_files[1]))
            elif key == 'TSV_File':
                items.append('file={}'.format(input_files[0]))
        if 'TSV_Oligos' in self._tool_inputs:
            items.append('oligos={}'.format(self._tool_inputs['TSV_Oligos'][0]))
        items.append('outputdir={}'.format(self._folder))
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the output file object with them
        :return: None
        """
        # As only one key is allowed, take the basename from the file(s) in the first key
        basename = super()._get_basename(list(self._tool_inputs.keys())[0], '.')
        self._tool_outputs['FASTA_Contig'] = [ToolIOFile(basename + '.trim.contigs.fasta')]
        self._tool_outputs['FASTA_Scrap'] = [ToolIOFile(basename + '.scrap.contigs.fasta')]
        self._tool_outputs['QUAL_Contig'] = [ToolIOFile(basename + '.trim.contigs.qual')]
        self._tool_outputs['QUAL_Scrap'] = [ToolIOFile(basename + '.scrap.contigs.qual')]
        self._tool_outputs['TSV_Report'] = [ToolIOFile(basename + '.contigs.report')]
        if 'TSV_File' in self._tool_inputs:
            self._tool_outputs['TSV_Groups'] = [ToolIOFile(basename + '.contigs.groups')]

    def _execute_tool(self):
        """
        Runs make.contigs
        :return: None
        """
        self._create_symlinks()
        self._build_command()
        self._execute_command()
        if self.__check_read_name_warning():
            self.__run_on_single_processor()
        self._symlink_cleanup()
        self._set_output()

    def __check_read_name_warning(self):
        """
        Function that checks whether a read name warning was given.
        :return: True if a read name warning was given
        """
        return bool(re.search('name mismatch in forward and reverse fastq file', self.stdout))

    def __run_on_single_processor(self):
        """
        Runs Mothur again but on a single processor regardless of the option specified in the database.
        Also removes the module load command as this is prepended in the run_tool method.
        :return: Nothing
        """
        self._command.command = re.sub(r'processors=(\d+)', 'processors=1', self._command.command)
        self._command.command = re.sub(r'module load .*;\s', '', self._command.command)
        self._execute_command()
