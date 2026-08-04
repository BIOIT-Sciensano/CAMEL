import re

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothur import Mothur


class MothurMakeContigs(Mothur):
    """
    The make.contigs command reads a forward fastq file and a reverse fastq file and outputs new fasta and report files.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_make_contigs')

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid:
        - Only one input key is allowed
        - Only FASTQ_PE, FASTA_PE, TSV_File and TSV_Oligos are allowed as keys
        - Only 2 files allowed for PE and 1 for TSV_File and TSV_Oligos files
        :return: None
        """
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidToolInputError(f'Too many input keys given for Mothur make.contigs: {self._tool_inputs}')
        self._required_input = list(self._tool_inputs.keys())
        super()._check_input()

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and input/output directories
        Example: ffastq=FileR1.fastq, rfastq=FileR2.fastq, inputdir=/test/data/input/,
        outputdir=/test/data/outputdir
        :return: String with the input parameters
        """
        items = []
        for key, input_files in self._tool_inputs.items():
            if key == 'FASTQ_PE':
                items.append(f'ffastq={input_files[0]}, rfastq={input_files[1]}')
            elif key == 'FASTA_PE':
                items.append(f'ffasta={input_files[0]}, rfasta={input_files[1]}')
            elif key == 'TSV_File':
                items.append(f'file={input_files[0]}')
        if 'TSV_Oligos' in self._tool_inputs:
            items.append(f"oligos={self._tool_inputs['TSV_Oligos'][0]}")
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the output file object with them
        :return: None
        """
        # As only one key is allowed, take the basename from the file(s) in the first key
        output_base = self._folder / next(iter(self._tool_inputs.values()))[0].basename
        self._tool_outputs['FASTA_Contig'] = [ToolIOFile(output_base.with_suffix('.trim.contigs.fasta'))]
        self._tool_outputs['FASTA_Scrap'] = [ToolIOFile(output_base.with_suffix('.scrap.contigs.fasta'))]
        self._tool_outputs['TSV_Report'] = [ToolIOFile(output_base.with_suffix('.contigs_report'))]
        if 'qfile' in self._parameters:
            self._tool_outputs['QUAL_Contig'] = [ToolIOFile(output_base.with_suffix('.trim.contigs.qual'))]
            self._tool_outputs['QUAL_Scrap'] = [ToolIOFile(output_base.with_suffix('.scrap.contigs.qual'))]
        if 'TSV_File' in self._tool_inputs:
            self._tool_outputs['TSV_Groups'] = [ToolIOFile(output_base.with_suffix('.contigs.groups'))]

    def _execute_tool(self) -> None:
        """
        Runs make.contigs
        :return: None
        """
        self._create_symlinks(self._temp_dir)
        self._build_command()
        self._execute_command()
        if self.__check_read_name_warning():
            self.__run_on_single_processor()
        self._set_output()

    def __check_read_name_warning(self) -> bool:
        """
        Function that checks whether a read name warning was given.
        :return: True if a read name warning was given
        """
        return bool(re.search('name mismatch in forward and reverse fastq file', self._command.stdout))

    def __run_on_single_processor(self) -> None:
        """
        Runs Mothur again but on a single processor regardless of the option specified in the database.
        Also removes the module load command as this is prepended in the run_tool method.
        :return: None
        """
        self._command.command = re.sub(r'processors=(\d+)', 'processors=1', self._command.command)
        self._command.command = re.sub(r'module load .*;\s', '', self._command.command)
        self._execute_command()
