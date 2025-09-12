import os.path

from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.qiime.qiime import Qiime


class QiimePickOpenReferenceOtus(Qiime):
    """
    Perform open-reference OTU picking
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('qiime_pick_open_reference_otus', '1.9.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA key is required
        - TSV_Taxonomy, FASTA_REF allowed as additional key
        - Only one file allowed per key
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('Invalid input files (keys) given for '
                                                 'pick_open_reference_otus: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.items():
            if key not in ['FASTA', 'FASTA_REF', 'TSV_Taxonomy']:
                raise InvalidToolInputError('Invalid input key given for pick_open_reference_otus: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidToolInputError('Invalid number (max = 1) of files in each key given for \
                                                     pick_open_reference_otus: {!r}'.format(self._tool_inputs))

    def _set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs['TSV_Otu'] = [ToolIOFile(self._folder / self.__get_otu_filename())]
        self._tool_outputs['BIOM'] = [ToolIOFile(self._folder / self.__get_otu_table_filename())]
        self._tool_outputs['FASTA_Ref'] = [ToolIOFile(self._folder / 'new_refseqs.fna')]
        self._tool_outputs['FASTA'] = [ToolIOFile(self._folder / 'rep_set.fna')]

    def __get_otu_filename(self):
        """
        Gets the name that the OTU map file will have
        :return: String with OTU map filename
        """
        otu_output = 'final_otu_map'
        if 'min_otu_size' in self._parameters:
            otu_output += '_mc{}.txt'.format(self._parameters['min_otu_size'].value)
        else:
            otu_output += '_mc2.txt'
        return otu_output

    def __get_otu_table_filename(self):
        """
        Gets the name that the OTU table file will have
        :return: String with OTU table filename
        """
        otu_table_output = 'otu_table'
        if 'min_otu_size' in self._parameters:
            otu_table_output += '_mc' + self._parameters['min_otu_size'].value
        else:
            otu_table_output += '_mc2'
        if 'suppress_taxonomy_assignment' not in self._parameters:
            otu_table_output += '_w_tax'
        if 'suppress_align_and_tree' not in self._parameters:
            otu_table_output += '_no_pynast_failures'
        otu_table_output += '.biom'
        return otu_table_output

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        input_string = '-i {}'.format(self._tool_inputs['FASTA'][0])
        if os.path.isfile(os.path.join(self._folder, 'parameters.txt')):
            input_string += ' -p {}'.format(os.path.join(self._folder, self._parameter_file))
        if 'FASTA_Ref' in self._tool_inputs:
            input_string += ' -r {}'.format(self._tool_inputs['FASTA_Ref'][0])
        input_string += ' -o {}'.format(self._folder)
        return input_string

    def _build_command(self):
        """
        Concatenates required parameters and options to build the command to run
        :return: None
        """
        options_string = self._build_options()
        options_string += ' --force'
        self._command.command = '{} {} {}'.format(self._tool_command, self._build_input_string(), options_string)
