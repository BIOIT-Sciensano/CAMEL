from camel.app.command.command import Command
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool
from camel.app.error import InvalidToolInputError, ToolExecutionError


class Vep(Tool):
    """
    ===========
    Variant Effect Predictor: annotates and predicts effects of variants in vcf.
    Oriented towards genetics more than somatic variants.
    ===========

    Required inputs:
    ----------------
    "VCF":              ToolIOFile object. VCF or TXT (whitespace seperated) file with variants.

    Output:
    -------
    "VCF":              ToolIOFile object. VCF file with annotations. Can also be txt file.
    "HTML":             ToolIOFile object. HTML file with summary of variants.

    Mandatory parameters:
    ---------------------
    - db_loc        Location of database to use for annotation
                    default: None
    - output_file   Output filename.
                    default:  vep_annotated.vcf

    Other parameters:
    -----------------
    - output_in_vcf    output in vcf format instead of default text format
    """

    def __init__(self) -> None:
        """
        Initialize Mutect1 tool.
        :return: None
        """
        super().__init__('Vep', '93')
        self._required_inputs = ['VCF']

    def _execute_tool(self) -> None:
        """
        Runs Mutect1
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self) -> None:
        """
        Check that input is valid (super method) and that required parameters are present.
        :return: None
        """
        for input_key in self._required_inputs:
            if input_key not in self._tool_inputs:
                raise InvalidToolInputError(self.name, f'Vep required {input_key} input is missing in tool inputs!')
        super()._check_input()

    def __build_command(self) -> None:
        """
        Build the command to run the tool.
        By default, use cache provided in DB_PATH and run offline.
        Concatenates strings into the command string.
        :return:
        """
        input_string = " -i {}".format(self._tool_inputs['VCF'][0].path)

        options_string = "--cache --offline "
        options_string += " ".join(self._build_options())

        self._command.command = ' '.join([self._tool_command, input_string, options_string])

    def __set_output(self) -> None:
        """
        Set the output specifications in the Camel ouptut list:
        - html file
        - vcf file
        :return: None
        """
        self._tool_outputs['OUT'] = [ToolIOFile(self._folder / self.get_param_value('output_file'))]

        html_file_name = self._folder / f"{self.get_param_value('output_file')}_summary.html"
        self._tool_outputs['HTML'] = [ToolIOFile(html_file_name)]

    def _check_command_output(self, command: Command) -> None:
        """
        Check the result of Vep tool run.
        :param command: Executed command
        :return: None
        """
        if not command.stdout == "" and "ERROR" in command.stdout:
            raise ToolExecutionError(self.name, f"Vep failed to run with message: \n{command.stdout}")
