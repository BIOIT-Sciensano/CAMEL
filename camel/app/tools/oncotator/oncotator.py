from camel.app.core.command import Command
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool
from camel.app.core.utils import toolutils


class Oncotator(Tool):
    """
    ===========
    Oncotator: anotates vcf with cancer-centered databases.
    ===========

    Required inputs:
    ----------------
    "IN_FILE":          ToolIOFile object. {TCGAMAF,SEG_FILE,VCF,MAFLITE}. See descriptions on https://gatkforums.broadinstitute.org/gatk/discussion/4177/what-input-files-can-i-annotate-with-oncotator

    Output:
    -------
    "OUT_FILE":         ToolIOFile object. {TCGAMAF,VCF,SIMPLE_TSV,TCGAVCF,SIMPLE_BED,GENE_LIST}

    Mandatory parameters:
    ---------------------
    - db_dir        Location of databases to use for annotation
                    default: None
    - output_file   Output filename.
                    default:  oncotator.tsv

    Other parameters:
    -----------------
    - canonical_tx_file:
                    File to use to supersede specific variants known to be clinically relevant
                    Default: None
    """

    def __init__(self) -> None:
        """
        Initialize Mutect1 tool.
        :return: None
        """
        super().__init__('Oncotator', '1.9.9.0')
        self._required_inputs = ['IN_FILE']

    def _execute_tool(self) -> None:
        """
        Runs tool
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_parameters(self) -> None:
        """
        Checks that parameters are valid.
        :return: None
        """
        super()._check_parameters()

    def _check_input(self) -> None:
        """
        Check that input is valid (super method) and that required parameters are present.
        :return: None
        """
        super()._check_input()

        for input_key in self._required_inputs:
            if input_key not in self._tool_inputs:
                raise InvalidToolInputError(f'Oncotator required {input_key} input is missing in tool inputs!')

    def __build_command(self) -> None:
        """
        Build the command to run the tool.
        :return:
        """
        input_string = f" {self._tool_inputs['IN_FILE'][0].path}"
        output_string = f" {self._parameters['output_file_name'].value}"
        options_string = " ".join(self._build_options(excluded_parameters=["output_file_name"]))
        self._command.command = f"{self._tool_command} {options_string} {input_string} {output_string} hg19"

    def __set_output(self) -> None:
        """
        Set the output specifications in the Camel output list.
        - output file
        :return: None
        """
        self._tool_outputs['OUT_FILE'] = [ToolIOFile(self._folder / self.get_param_value('output_file_name'))]

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
