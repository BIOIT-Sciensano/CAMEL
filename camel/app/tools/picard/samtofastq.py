from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.picard.picard import Picard


class SamToFastq(Picard):
    """
    Class for Picard FastqToSam function
    """

    def __init__(self):
        """
        Initialize a picard tool
        :return: None
        """
        super().__init__('Picard SamToFastq', '2.23.3')
        self._output_type = 'FASTQ'
        self._specific_parameters = ['paired_end', 'output']

    def _set_output(self) -> None:
        """
        Set the output specification, depending on paired-end parameter
        This function overrules the function in the parent class.
        :return: None
        """
        self.__set_output_string()

        if self._parameters['paired_end']:
            self._tool_outputs[self._output_type] = [
                ToolIOFile(Path(self._folder) / f"{self._parameters['output'].value}_R1.fastq"),
                ToolIOFile(Path(self._folder) / f"{self._parameters['output'].value}_R2.fastq")
            ]
        else:
            self._tool_outputs[self._output_type] = [
                ToolIOFile(Path(self._folder) / f"{self._parameters['output'].value}_R1.fastq")
            ]

    def __set_output_string(self) -> None:
        """
        Sets the output string. If paired-end, --FASTQ and --SECOND_END_FASTQ need to be defined
        :return: None
        """
        if self._parameters['paired_end']:
            self._output_string = f"FASTQ={self._parameters['output'].value}_R1.fastq " \
                                  f"SECOND_END_FASTQ={self._parameters['output'].value}_R2.fastq"
        else:
            self._output_string = f"FASTQ={self._parameters['output'].value}_R1.fastq "
