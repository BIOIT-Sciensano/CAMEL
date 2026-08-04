from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.seqtk.seqtk import Seqtk


class SeqtkConvert(Seqtk):
    """
    Class that converts fastq into fasta file using seqtk
    """

    def __init__(self) -> None:
        """
        Initialize seqtk convert
        :return: None
        """
        super().__init__('Seqtk Convert', '1.4')
        self._function_name = 'Convert'
        self._supported_inputs = ['FASTQ']
        self._specific_parameters = ['output_file']

    def _get_input_string(self) -> Path:
        """
        Returns the input specification
        :return: input_string containing input specification
        """
        return self._tool_inputs['FASTQ'][0].path

    def _set_output(self) -> None:
        """
        Set the output specification
        :return: None
        """
        self._output_string = self._folder / self._parameters['output_file'].value
        self._tool_outputs['FASTA'] = [ToolIOFile(self._output_string)]
