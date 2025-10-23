import os

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.mothur.mothur import Mothur


class MothurRemoveSeqs(Mothur):
    """
    The remove.seqs command takes a list of sequence names and either a fastq, fasta, name, group, list, count or
    align.report file to generate a new file that does not contain the sequences in the list.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_remove_seqs', '1.39.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - TSV_Accnos is required
        - In addition only two extra keys are permitted
        - Possible additional key comes from: 'FASTA', 'TSV_Names', 'TSV_Counts', 'TSV_Groups',
          'TSV_AlignReport', 'TSV_List', 'TSV_Taxonomy', 'TSV_Qfile', 'FASTQ'
        - Only one file per key is allowed
        :return: None
        """
        super()._check_input()
        if 'TSV_Accnos' not in self._tool_inputs:
            raise InvalidToolInputError('Not enough valid input files given for Mothur '
                                                 'remove.seqs: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) > 3:
            raise InvalidToolInputError('Too many input keys given for Mothur remove.seqs: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.items():
            if key not in ['TSV_Accnos', 'FASTA', 'TSV_Names', 'TSV_Counts', 'TSV_Groups',
                           'TSV_AlignReport', 'TSV_List', 'TSV_Taxonomy', 'TSV_Qfile', 'FASTQ']:
                raise InvalidToolInputError('Invalid input key given for Mothur remove.seqs: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidToolInputError('Invalid number (max = 1) of files given for Mothur \
                                                     remove.seqs: {!r}'.format(self._tool_inputs))
        self.__check_empty_input()

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        # TSV_Accnos is required so will always be available
        items = ['accnos={}'.format(self._tool_inputs['TSV_Accnos'][0])]
        input_parameters = {'FASTA': 'fasta=',
                            'FASTQ': 'fastq=',
                            'TSV_Names': 'name=',
                            'TSV_Counts': 'count=',
                            'TSV_Groups': 'group=',
                            'TSV_AlignReport': 'alignreport=',
                            'TSV_List': 'list=',
                            'TSV_Taxonomy': 'taxonomy=',
                            'TSV_Qfile': 'qfile='}
        for key, input_files in self._tool_inputs.items():
            # Only two keys are possible so the one that is not TSV_Accnos
            # will define the option flag that is needed
            if key != 'TSV_Accnos':
                items.append('{}{}'.format(input_parameters[key], input_files[0].path))
        items.append('outputdir={}'.format(self._folder))
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        # The extension of the output files is not always built in the same way. A dictionary
        # holds the suffixes that can be used when getting the base name
        output_extensions = {'FASTA': ['.', '.pick.fasta'],
                             'FASTQ': ['.', '.pick.fastq'],
                             'TSV_Names': ['.', '.pick.names'],
                             'TSV_Counts': ['.', '.pick.count_table'],
                             'TSV_Groups': ['.', '.pick.groups'],
                             'TSV_AlignReport': ['.align.report', '.pick.align.report'],
                             'TSV_List': ['.', '.pick.list'],
                             'TSV_Taxonomy': ['.', '.pick.taxonomy'],
                             'TSV_Qfile': ['.', '.pick.qual']}
        for key, input_files in self._tool_inputs.items():
            # The TSV_Accnos file does not directly lead to output
            if key != 'TSV_Accnos':
                basename = super()._get_basename(key, output_extensions[key][0])
                self._tool_outputs[key] = [ToolIOFile(basename + output_extensions[key][1])]

    def __check_empty_input(self):
        if os.path.getsize(self._tool_inputs['TSV_Accnos'][0].path) == 0:
            with open(self._tool_inputs['TSV_Accnos'][0].path, 'wt', encoding='utf-8') as outf:
                outf.write('adding_dummy_entry_as_original_file_is_empty\n')
                logger.warning('WARNING: ACCNOS file was empty, added a dummy record!')
