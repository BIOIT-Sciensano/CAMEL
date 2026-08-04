import re
from pathlib import Path
from typing import Optional

from camel.app.core.errors import InvalidParameterError, InvalidToolInputError
from camel.app.tools.picard.picard import Picard


class FastqToSam(Picard):
    """
    Class for Picard FastqToSam function
    """

    def __init__(self):
        """
        Initialize a picard tool
        :return: None
        """
        super().__init__('Picard FastqToSam', '2.23.3')

        self._required_inputs = []

    def _set_input(self) -> None:
        """
        Set input for the picard function
        :return: None
        """
        if 'FASTQ_PE' in self._tool_inputs:
            self._input_string = f"FASTQ={self._tool_inputs['FASTQ_PE'][0].path} " \
                                 f"FASTQ2={self._tool_inputs['FASTQ_PE'][1].path}"
        elif 'FASTQ_SE' in self._tool_inputs:
            self._input_string = f"FASTQ={self._tool_inputs['FASTQ_SE'][0].path}"
        else:
            InvalidToolInputError('Picard FastqToSam requires FASTQ_SE or FASTQ_PE input.')

        if 'SAMPLE_NAME' in self._tool_inputs:
            # if SAMPLE_NAME specified, it will replace the default values of parameters: RG_sample_name, RG_name
            self._specific_parameters = ['RG_sample_name', 'RG_name']
            self._input_string += f" SM={self._tool_inputs['SAMPLE_NAME'][0].value} RG={self._tool_inputs['SAMPLE_NAME'][0].value}"

    def _set_output(self) -> None:
        """
        Set output for FastqToSam. Extension determines whether a SAM or BAM file will be made
        :return: None
        """
        self._output_type = Path(self._parameters['output'].value).suffix.strip(".").upper()
        if self._output_type not in ["SAM", "BAM"]:
            raise InvalidParameterError("Picard FastqToSam: output file extension should be .bam or .sam")
        super()._set_output()

    def _set_informs(self, stderr: Optional[str] = None) -> None:
        """
        Analyse the result of picard run and update tool.informs.
            - Quality encoding: Quality score format (depending on sequencing technology)
            - Reads processed: Number of reads processed
        :return: None
        """
        for line in (self._command.stderr if stderr is None else stderr).splitlines():
            m = re.search('Auto-detected quality format as: ([a-zA-Z]+)', line)
            if m:
                if m.group(1) == 'Standard':
                    self.informs['quality_encoding'] = 'phred33'
                elif m.group(1) == 'Illumina':
                    self.informs['quality_encoding'] = 'phred64'
                elif m.group(1) == 'Solexa':
                    self.informs['quality_encoding'] = 'solexa66'
                else:
                    self.informs['quality_encoding'] = f'[UNKNOWN]{m.group(1)}'

            m = re.search(r'Processed (\d+) fastq reads', line)
            if m:
                self.informs['reads_processed'] = m.group(1)
