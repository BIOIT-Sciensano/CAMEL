from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.bedtools.bedtools import Bedtools


class BedtoolsMerge(Bedtools):
    """
    Tool class for Bedtools merge function.
    ==========================
    Bedtools merge 2.25.0
    ==========================
    https://bedtools.readthedocs.io/en/latest/content/tools/merge.html?highlight=merge
    bedtools merge combines overlapping or “book-ended” features in an interval file into a single feature which spans
    all of the combined features.
    Can use BAM as input and generate BED as output.

    This camel implementation will generate BED file if BAM is given as input.

    Required inputs:
    ----------------
    'BAM'/'BED':        Input BAM or BED file. (Max one file at a time)

    Output:
    -------
    'BED':              Bed file with regions covered by input BAM/BED file.

    Mandatory parameters:
    ---------------------
    - output_filename   Default value: 'output.bed'
    """

    def __init__(self) -> None:
        """
        Initialize a bedtools tool.
        :return: None
        """
        super().__init__('bedtools merge', '2.31.0')
        self.__input_type = ""

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        if self.__input_type == "BAM":
            self.update_parameters(use_bed_output=True)

        self.__set_output()
        self.__build_command()
        self._execute_command()

    def __build_command(self):
        """
        Builds the command with input, options and output strings.
        :return: None
        """
        input_string = f"-i {self._tool_inputs[self.__input_type][0].path}"
        build_options = ' '.join(self._build_options(excluded_parameters=['output_filename']))
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
        super()._check_input()

    def _check_required_inputs(self):
        """
        Checks that required input is present and sets the input type.
        Supersedes the same function in the Bedtools class.
        Inputs must be one of either BAM or BED.
        :return: None
        """
        if len(self._tool_inputs) != 1:
            raise InvalidToolInputError(
                f"{len(self._tool_inputs)} input file(s) specified. Bedtools merge takes exactly ONE input (BAM or BED).")
        elif "BAM" in self._tool_inputs:
            self.__input_type = "BAM"
        elif "BED" in self._tool_inputs:
            self.__input_type = "BED"
        else:
            raise InvalidToolInputError(
                f"Input file specified with wrong file type ({list(self._tool_inputs.keys())[0]}). Accepted types are BAM or BED.")

    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['BED'] = [ToolIOFile(self._folder / self._parameters['output_filename'].value)]
