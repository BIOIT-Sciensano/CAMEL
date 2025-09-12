from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.seqtk.seqtk import Seqtk


class SeqtkMergePE(Seqtk):
    """
    Interleaves two PE FASTA/Q files.
    """

    def __init__(self) -> None:
        """
        Initialize seqtk convert
        :return: None
        """
        super().__init__('Seqtk mergepe', '1.4')
        self._function_name = 'mergepe'
        self._supported_inputs = ['FASTQ_PE']
        self._specific_parameters = ['output_file']

    def _get_input_string(self) -> str:
        """
        Returns the input specification
        :return: input_string containing input specification
        """
        return ' '.join(str(x) for x in self._tool_inputs['FASTQ_PE'])

    def _set_output(self) -> None:
        """
        Set the output specification
        :return: None
        """
        self._output_string = self._folder / self._parameters['output_file'].value
        self._tool_outputs['FASTQ'] = [ToolIOFile(self._output_string)]
