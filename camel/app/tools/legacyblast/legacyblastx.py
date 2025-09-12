import os
from pathlib import Path

from Bio import SeqIO

from camel.app.command.command import Command
from camel.app.components.files.fileutils import FileUtils
from camel.app.error import InvalidToolInputError, ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class LegacyBlastx(Tool):
    """
    Map a nucleotide sequence against a protein sequence database
    """

    def __init__(self):
        """
        Initialize tool
                :return: None
        """
        super().__init__('legacy_blastx', '2.2.22')
        self._fasta = None

    def _execute_tool(self):
        """
        Runs Carma
        :return: None
        """
        self.__prepare_input_files()
        for infile in self._fasta:
            self.__build_command(infile)
            self._execute_command()
        self.__set_output()
        if len(self._fasta) > 1:
            self.__concatenate_split_outputs()

    def __concatenate_split_outputs(self):
        """
        Concatenates the outputs of the split input file into one big output file
        :return: None
        """
        FileUtils.concatenate_files(self._tool_outputs['BLASTX'][0].path,
                                    [os.path.splitext(f.path)[0] + '.blastx' for f in self._fasta])
        for infile in self._fasta:
            os.remove(infile.path)
            os.remove(os.path.splitext(infile.path)[0] + '.blastx')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA and DB key are required
        - Only one input file allowed
        - No other input keys are allowed
        :return: None
        """
        super(LegacyBlastx, self)._check_input()
        if 'FASTA' not in self._tool_inputs or 'DB' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA or DB input files missing from input for '
                                                 'legacy blastx: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs['FASTA']) != 1 or len(self._tool_inputs['DB']) != 1:
            raise InvalidToolInputError('Invalid number (max = 1) of files per key given '
                                                 'for legacy blastx: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) > 2:
            raise InvalidToolInputError('Too many input keys given for legacy blastx '
                                                 '(only FASTA and DB allowed): {!r}'.format(self._tool_inputs))

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        if 'ERROR' in command.stderr:
            raise ToolExecutionError(self.name, "Command execution failed (stderr: {}).".format(command.stderr))
        if self._command.exit_code != 0:
            raise ToolExecutionError(self.name, "Command execution failed (Exit code: {})".format(command.exit_code))

    def __get_basename(self):
        """
        Returns the prefix that will be used in the output.
        :return: String with the prefix used in the output
        """
        infile = os.path.basename(self._tool_inputs['FASTA'][0].path)
        return os.path.join(self._folder, os.path.splitext(infile)[0])

    def __set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        basename = self.__get_basename()
        self._tool_outputs['BLASTX'] = [ToolIOFile(Path(basename + '.blastx'))]

    def __build_input_string(self, infile):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = ['-p blastx',
                 '-d {}'.format(self._tool_inputs['DB'][0].path),
                 '-i {}'.format(infile.path),
                 '-o {}.blastx'.format(os.path.join(self._folder, os.path.splitext(os.path.basename(infile.path))[0]))]
        return ' '.join(items)

    def __prepare_input_files(self):
        """
        Prepares the input files by splitting a large fasta file into smaller chunks if necessary
        :return:
        """
        if self.__is_large_fasta():
            self.__split_fasta()
        else:
            self._fasta = [self._tool_inputs['FASTA'][0]]

    def __build_command(self, infile):
        """
        Concatenates required parameters and options to build the command
        :param infile: The name of the input file
        :return: None
        """
        input_string = self.__build_input_string(infile)
        options_string = ' '.join(self._build_options())
        self._command.command = '{} {} {}'.format(self._tool_command, input_string, options_string)

    def __is_large_fasta(self):
        """
        Checks whether the input FASTA file has more than 500 records
        :return: True when there are more than 500 records in the FASTA file
        """
        return len(SeqIO.index(self._tool_inputs['FASTA'][0].path, 'fasta')) > 500

    @staticmethod
    def __create_iterator_batch(iterator, batch_size):
        """
        Returns lists of length batch_size.

        This can be used on any iterator, for example to batch up
        SeqRecord objects from Bio.SeqIO.parse(...), or to batch
        Alignment objects from Bio.AlignIO.parse(...), or simply
        lines from a file handle.

        This is a generator function, and it returns lists of the
        entries from the supplied iterator.  Each list will have
        batch_size entries, although the final list may be shorter.
        """
        entry = True
        while entry:
            batch = []
            while len(batch) < batch_size:
                try:
                    entry = next(iterator)
                except StopIteration:
                    entry = None
                if entry is None:
                    break
                batch.append(entry)
            if batch:
                yield batch

    def __split_fasta(self):
        """
        Splits the input FASTA file into smaller chuncks in order to speed up the blastx analysis
        :return:
        """
        self._fasta = []
        record_iter = SeqIO.parse(self._tool_inputs['FASTA'][0].path, 'fasta')
        for i, batch in enumerate(self.__create_iterator_batch(record_iter, 100)):
            with open(os.path.join(self._folder, 'temp_{}.fasta'.format(i)), 'wb') as f:
                SeqIO.write(batch, f, 'fasta')
                self._fasta.append(ToolIOFile(self._folder / f'temp_{i}.fasta'))
