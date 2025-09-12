from camel.app.error import InvalidToolInputError
from camel.app.tools.bedtools.bedtools import Bedtools
from camel.app.io.tooliofile import ToolIOFile


class BedtoolsBamToBed(Bedtools):
    """
    Tool class for Bedtools BamtoBed function.

    ==========================
    Bedtools bamtobed 2.25.0
    ==========================
    https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html
    bedtools bamtobed is a conversion utility that converts sequence alignments in BAM format into BED, BED12, and/or
    BEDPE records.
    This camel implementation only allows conversion to BED format.

    Required inputs:
    ----------------
    'BAM':              Input bam file. (Max one file at a time)

    Output:
    -------
    'BED':              Bed file with regions covered by input BAM file.

    Mandatory parameters:
    ---------------------
    - output_filename   Default value: 'output.bed'
    """

    def __init__(self):
        """
        Initialize a bedtools tool.
        :return: None
        """
        super().__init__('bedtools bamtobed', '2.31.0')
        self._required_inputs = ['BAM']

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __build_command(self):
        """
        Builds the command with input, options and output strings.
        :return: None
        """
        build_options = ' '.join(self._build_options(excluded_parameters=['output_filename']))
        input_string = "-i {} ".format(self._tool_inputs['BAM'][0].path)
        output_string = '> ' + self._parameters['output_filename'].value
        self._command.command = ' '.join([
            self._tool_command,
            build_options,
            input_string,
            output_string])

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        self._check_required_inputs()
        if len(self._tool_inputs['BAM']) != 1:
            raise InvalidToolInputError("Exactly one BAM input file expected.")
        super()._check_input()

    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['BED'] = [ToolIOFile(self._folder / self._parameters['output_filename'].value)]
