from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.qiime.qiime import Qiime


class QiimeJoinPairedEnds(Qiime):
    """
    This script takes forward and reverse Illumina reads and joins them using the method chosen.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('qiime_join_paired_ends', '1.9.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTQ_PE key is required
        - No additional keys allowed
        - Only two files allowed per key
        :return: None
        """
        if 'FASTQ_PE' not in self._tool_inputs:
            raise InvalidToolInputError('Invalid input files (keys) given for join_paired_ends: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidToolInputError('Too many input keys given for join_paired_ends: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs['FASTQ_PE']) != 2:
            raise InvalidToolInputError('Invalid number (!= 2) of files in each key given for \
                                                  join_paired_ends: {!r}'.format(self._tool_inputs))

    def _set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs['FASTQ'] = [ToolIOFile(self._folder / 'fastqjoin.join.fastq')]
        self._tool_outputs['FASTQ_Unjoined'] = [
            ToolIOFile(self._folder / 'fastqjoin.un1.fastq'),
            ToolIOFile(self._folder / 'fastqjoin.un2.fastq')
        ]

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        return ' '.join([
            f'-f {self._tool_inputs["FASTQ_PE"][0]}',
            f'-r {self._tool_inputs["FASTQ_PE"][1]}',
            f'-o {self._folder}'
        ])
