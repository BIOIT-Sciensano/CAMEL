import os.path

from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.qiime.qiime import Qiime


class QiimePickClosedReferenceOtus(Qiime):
    """
    Perform open-reference OTU picking
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('qiime_pick_closed_reference_otus', '1.9.1')

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
                                                 'pick_closed_reference_otus: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.items():
            if key not in ['FASTA', 'FASTA_Ref', 'TSV_Taxonomy']:
                raise InvalidToolInputError('Invalid input key given for '
                                                     'pick_closed_reference_otus: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidToolInputError('Invalid number (max = 1) of files in each key given for \
                                                     pick_closed_reference_otus: {!r}'.format(self._tool_inputs))

    def _set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs['BIOM'] = [ToolIOFile(self._folder / 'otu_table.biom')]

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
        if 'TSV_Taxonomy' in self._tool_inputs:
            input_string += ' -t {}'.format(self._tool_inputs['TSV_Taxonomy'][0])
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
