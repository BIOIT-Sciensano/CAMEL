from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurSubSample(Mothur):
    """
    The sub.sample command can be used as a way to normalize your data, or create a smaller set from your original set.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_sub_sample', version=None)

    def _check_input(self):
        """
        Checks whether the given inputs are valid
        Some keys are considered primary input and cannot be combined with other primary input keys.
        :return: None
        """
        super()._check_input()
        allowed_primary_input = ['FASTA', 'TSV_List', 'TSV_Shared', 'TSV_Rabund', 'TSV_Sabund']
        for key, input_files in self._tool_inputs.items():
            if key not in allowed_primary_input and key not in ['TSV_Groups', 'TSV_Counts', 'TSV_Names']:
                raise InvalidToolInputError('Invalid input key given for Mothur sub.sample: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidToolInputError('Invalid number (max = 1) of files given for Mothur \
                                                     sub.sample: {!r}'.format(self._tool_inputs))
        # Check if more than one of the primary input keys is present
        if len(set(allowed_primary_input).intersection(self._tool_inputs.keys())) > 1:
            raise InvalidToolInputError('Too many primary input keys given for '
                                                 'Mothur sub.sample: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = []
        if 'FASTA' in self._tool_inputs:
            items.append('fasta={}'.format(self._tool_inputs['FASTA'][0]))
        elif 'TSV_List' in self._tool_inputs:
            items.append('list={}'.format(self._tool_inputs['TSV_List'][0]))
        elif 'TSV_Shared' in self._tool_inputs:
            items.append('shared={}'.format(self._tool_inputs['TSV_Shared'][0]))
        elif 'TSV_Rabund' in self._tool_inputs:
            items.append('rabund={}'.format(self._tool_inputs['TSV_Rabund'][0]))
        elif 'TSV_Sabund' in self._tool_inputs:
            items.append('sabund={}'.format(self._tool_inputs['TSV_Sabund'][0]))
        if 'TSV_Names' in self._tool_inputs:
            items.append('name={}'.format(self._tool_inputs['TSV_Names'][0]))
        if 'TSV_Counts' in self._tool_inputs:
            items.append('count={}'.format(self._tool_inputs['TSV_Counts'][0]))
        if 'TSV_Groups' in self._tool_inputs:
            items.append('group={}'.format(self._tool_inputs['TSV_Groups'][0]))
        items.append('outputdir={}'.format(self._folder))
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them.
        :return: None
        """
        output_extensions = {'FASTA': ['.', '.subsample.fasta'],
                             'TSV_List': ['.', '.subsample.list'],
                             'TSV_Shared': ['.', '.subsample.shared'],
                             'TSV_Rabund': ['.', '.subsample.rabund'],
                             'TSV_Sabund': ['.', '.subsample.sabund']}
        for key, input_files in self._tool_inputs.items():
            if key in output_extensions:
                basename = super()._get_basename(key, output_extensions[key][0])
                self._tool_outputs[key] = []
                if self.__get_labels() is not None:
                    for label in self.__get_labels():
                        self._tool_outputs[key] += [ToolIOFile(basename + '.' + label + output_extensions[key][1])]
                else:
                    self._tool_outputs[key] += [ToolIOFile(basename + output_extensions[key][1])]

    def __get_labels(self):
        """
        Returns the labels that are specified as a parameter.
        :return: List of labels
        """
        if 'label' in self._parameters:
            return self._parameters['label'].value.strip().split('-')
        return None
