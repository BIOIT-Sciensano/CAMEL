from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.qiime.qiime import Qiime


class QiimeSplitLibrariesFastq(Qiime):
    """
    Split_libraries_fastq demultiplexes fastq files if needed and performs several quality controls
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('qiime_split_libraries_fastq', '1.9.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA or FASTQ key is required
        - Only TSV_Map allowed as additional key
        - Only one file allowed per key
        :return: None
        """
        if 'FASTA' in self._tool_inputs == 'FASTQ' in self._tool_inputs:
            raise InvalidToolInputError('Invalid input files (keys) given for '
                                                 f'split_libraries_fastq: {self._tool_inputs!r}')
        for key, input_files in self._tool_inputs.items():
            if key not in ['FASTA', 'FASTQ', 'TSV_Map']:
                raise InvalidToolInputError(f'Invalid input key given for split_libaries_fastq: {self._tool_inputs!r}')
            if len(input_files) != 1:
                raise InvalidToolInputError(f'Invalid number (max = 1) of files in each key given for \
                                                     split_libraries_fastq: {self._tool_inputs!r}')

    def _set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs['TSV_Histogram'] = [ToolIOFile(self._folder / 'histograms.txt')]
        self._tool_outputs['FASTA'] = [ToolIOFile(self._folder / 'seqs.fna')]
        self._tool_outputs['LOG'] = [ToolIOFile(self._folder / 'split_library_log.txt')]

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        if 'FASTA' in self._tool_inputs:
            input_string = '-i {}'.format(self._tool_inputs['FASTA'][0])
        else:
            input_string = '-i {}'.format(self._tool_inputs['FASTQ'][0])
        input_string += f' -o {self._folder}'
        return input_string
