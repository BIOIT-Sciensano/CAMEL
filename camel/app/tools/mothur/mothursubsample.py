from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothur import Mothur


class MothurSubSample(Mothur):
    """
    The sub.sample command can be used as a way to normalize your data, or create a smaller set from your original set.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_sub_sample')
        self._optional_input = ['TSV_Groups', 'TSV_Counts', 'TSV_Names']

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid
        Some keys are considered primary input and cannot be combined with other primary input keys.
        :return: None
        """
        allowed_primary_input = {'FASTA', 'TSV_List', 'TSV_Shared', 'TSV_Rabund', 'TSV_Sabund'}
        # Check if more than one of the primary input keys is present
        if len(set(allowed_primary_input).intersection(self._tool_inputs.keys())) > 1:
            raise InvalidToolInputError(f'Too many primary input keys given for Mothur sub.sample: {self._tool_inputs}')
        else:
            self._required_input = list(set(allowed_primary_input).intersection(self._tool_inputs.keys()))
        super()._check_input()

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = []
        if 'FASTA' in self._tool_inputs:
            items.append(f"fasta={self._tool_inputs['FASTA'][0]}")
        elif 'TSV_List' in self._tool_inputs:
            items.append(f"list={self._tool_inputs['TSV_List'][0]}")
        elif 'TSV_Shared' in self._tool_inputs:
            items.append(f"shared={self._tool_inputs['TSV_Shared'][0]}")
        elif 'TSV_Rabund' in self._tool_inputs:
            items.append(f"rabund={self._tool_inputs['TSV_Rabund'][0]}")
        elif 'TSV_Sabund' in self._tool_inputs:
            items.append(f"sabund={self._tool_inputs['TSV_Sabund'][0]}")
        if 'TSV_Names' in self._tool_inputs:
            items.append(f"name={self._tool_inputs['TSV_Names'][0]}")
        if 'TSV_Counts' in self._tool_inputs:
            items.append(f"count={self._tool_inputs['TSV_Counts'][0]}")
        if 'TSV_Groups' in self._tool_inputs:
            items.append(f"group={self._tool_inputs['TSV_Groups'][0]}")
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them.
        :return: None
        """
        output_extensions = {'FASTA': '.subsample.fasta',
                             'TSV_List': '.subsample.list',
                             'TSV_Shared': '.subsample.shared',
                             'TSV_Rabund': '.subsample.rabund',
                             'TSV_Sabund': '.subsample.sabund'}
        for key in self._tool_inputs.keys():
            if key in output_extensions:
                basename = self._get_basename(key)
                self._tool_outputs[key] = []
                if self.__get_labels() is not None:
                    for label in self.__get_labels():
                        self._tool_outputs[key] += [ToolIOFile(basename.with_suffix(f'.{label}{output_extensions[key]}'))]
                else:
                    self._tool_outputs[key] += [ToolIOFile(basename.with_suffix(output_extensions[key]))]

    def __get_labels(self) -> list[str] | None:
        """
        Returns the labels that are specified as a parameter.
        :return: List of labels
        """
        if 'label' in self._parameters:
            return self.get_param_value('label').strip().split('-')
        return None
