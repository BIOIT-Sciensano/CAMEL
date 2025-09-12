from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.star.star import Star


class StarAlign(Star):
    """
    Align spliced transcripts with STAR.
    """

    def __init__(self) -> None:
        """
        Initializes STAR 2.7.11b alignment.
        :return: None
        """
        super().__init__('STAR', '2.7.11b')
        self._required_inputs = ['FASTQ', 'INDEX_DIR']
        self._input_string = "--runMode alignReads"

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        super()._check_input()
        if len(self._tool_inputs['FASTQ']) == 0:
            raise InvalidToolInputError("No FASTQ input provided")
        elif len(self._tool_inputs['FASTQ']) > 2:
            raise InvalidToolInputError("Too many inputs of type FASTQ")

    def _set_input(self) -> None:
        """
        Sets the input specification and the input string.
        :return: None
        """
        option_fastq = "--readFilesIn"
        for fastq in self._tool_inputs['FASTQ']:
            option_fastq += f" {fastq.path}"

        option_gtf = ""
        if 'GTF' in self._tool_inputs:
            option_gtf = f" --sjdbGGTFfile {self._tool_inputs['GTF'][0].path}"

        option_index_dir = f"--genomeDir {self._tool_inputs['INDEX_DIR'][0].path}"

        self._input_string = " ".join([self._input_string,
                                       option_fastq,
                                       option_gtf,
                                       option_index_dir])

    def _set_output(self) -> None:
        """
        Sets the output specification and the output string.
        :return: None
        """
        filename_output = self._parameters['filename_output'].value
        self._output_string += f" > {filename_output}"
        self._tool_outputs['SAM'] = [ToolIOFile(self.folder / filename_output)]
