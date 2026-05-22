import re

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.tool import Tool


class Mutect1(Tool):
    """
    ===========
    Mutect (v1).
    ===========
    Performs somatic variant calling, e.g. for oncology-related NGS data. 
    Mutect V1 only calls snps, not indels. For indels, use Mutect2.
    
    Required inputs:
    ----------------
    "BAM_TUMOR":        ToolIOFile object. BAM file with tumour data
    
    Optional input:
    ---------------
    "BAM_NORMAL":       ToolIOFile object. BAM file with normal data for tumor-normal matching.
    "FASTA_REF":        ToolIOFile object. FASTA file containing the reference genome.
    "VCF_DBSNP":        ToolIOFile object. DbSNP reference vcf file location.
    "TXT_intervals":    ToolIOFile object. Intervals list to restrict search by GATK. Accelerates analysis. Bed or GATK intervals list 
    
    Output:
    -------
    "TXT_CALL_STATS":   ToolIOFile object. GATK Call stats text based file. Parse-able by scripts or in excel sheets.
    "VCF":              ToolIOFile object. VCF file.
    
    Mandatory parameters:
    ---------------------
    - output_callstats_file
                    default value:  call_stats.txt
    """

    def __init__(self):
        """
        Initialize Mutect1 tool.
        :return: None
        """
        super().__init__('Mutect1', '1.1.7')
        self._required_inputs = ['BAM_TUMOR', 'FASTA_REF']

    def _execute_tool(self):
        """
        Runs Mutect1
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Check that input is valid (super method) and that required parameters are present.
        :return: None
        """
        super(Mutect1, self)._check_input()

        for input_key in self._required_inputs:
            if input_key not in self._tool_inputs:
                raise InvalidToolInputError(
                    f'Mutect1 required {input_key} input is missing in tool inputs!')

    def _check_parameters(self):
        """
        Checks that parameters are valid.
        :return: None 
        """
        super(Mutect1, self)._check_parameters()

    def __build_command(self):
        """
        Build the command to run the tool.
        Concatenates strings into the command string.
        :return: 
        """
        input_string = self.__create_input_string()
        options_string = ' '.join(self._build_options())
        self._command.command = ' '.join([self._tool_command, input_string, options_string])

    def __create_input_string(self):
        """
        Add the input specification in the input_string: 
        - DBSNP VCF file
        - reference fasta (default or superseded)
        - tumour BAM file
        - normal tissue BAM file (optional)
        - Intervals list for acceleration
        :return: Input_string
        """
        input_string = ""
        # set reference genome
        input_string += "-R {} ".format(self._tool_inputs['FASTA_REF'][0].path)

        # Use intervals to restrict search if supplied.
        if 'TXT_intervals' in self._tool_inputs:
            input_string += "-L {} ".format(self._tool_inputs['TXT_intervals'][0].path)

        # set reference dbSNP db
        if 'VCF_DBSNP' in self._tool_inputs:
            input_string += "--dbsnp {} ".format(self._tool_inputs['VCF_DBSNP'][0].path)

        input_string += "-I:tumor {} ".format(self._tool_inputs['BAM_TUMOR'][0].path)

        if 'BAM_NORMAL' in self._tool_inputs:
            input_string += "-I:normal {} ".format(self._tool_inputs['BAM_NORMAL'][0].path)

        return input_string

    def __set_output(self):
        """
        Set the output specifications in the Camel ouptut list: 
        - call_stats file
        - vcf file
        :return: None
        """
        self._tool_outputs['TXT_CALL_STATS'] = [ToolIOFile(self._folder / self.get_param_value('output_callstats_file'))]
        self._tool_outputs['VCF'] = [ToolIOFile(self._folder / self.get_param_value('output_vcf_file'))]

    def _check_command_output(self, command: Command) -> None:
        """
        Check the result of Mutect1 tool run.
        :param command: Command
        :return: None
        """
        if len(self._command.stdout.split('\n')) > 1:
            if not re.match('Exit status: 0', self.stdout.split('\n')[-2].rstrip()):
                raise ToolExecutionError(f"Mutect1 fails to run, message: \n{self.stdout}")
