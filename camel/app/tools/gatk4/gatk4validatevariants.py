import shutil
from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.errors import ToolExecutionError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.gatk4.gatk4 import GATK4


class GATK4ValidateVariants(GATK4):
    """
    Class for GATK4 ValidateVariants function, validate the adherence of a file to VCF format.

    Required inputs:
    ----------------
    'VCF':              ToolIOFile object. Input VCF file.
    'FASTA_REF':        ToolIOFile object. FASTA file containing the reference genome.

    """

    def __init__(self) -> None:
        """
        Initialize GATK4 ValidateVariants
        :return: None
        """
        super().__init__('gatk4 ValidateVariants', '4.1.9.0')

        self._required_inputs = ["VCF", "FASTA_REF"]
        self._specific_parameters = ["output"]

    def _set_input(self) -> None:
        """
        Set the input specification. Overrules method in parent class
        :return: None
        """
        super()._set_input()
        if 'VCF_dbsnp' in self._tool_inputs:
            self._input_string += f"--dbsnp {self._tool_inputs['VCF_dbsnp'][0].path} "

    def _set_output(self) -> None:
        """
        Set the output specification. Overrules method in parent class
        :return: None
        """
        self._tool_outputs['TXT_metrics'] = [ToolIOFile(Path(self._folder) / self._parameters['output'].value)]

    def _print_output(self) -> None:
        """
        This function does not create an output file, but outputs everything to stderr
        Print output file for snakemake's sake
        :return: None
        """
        with open(self._tool_outputs['TXT_metrics'][0].path, 'w') as output_file:
            output_file.write(self._command.stderr)

    def _check_command_output(self, command: Command) -> None:
        """
        Check the result. Error in the VCF file should not automatically cause the pipeline to terminate.
        :return: None
        """
        try:
            super()._check_command_output(command)
        except ToolExecutionError:
            logger.warning(f"GATK tool {self._name} failed to run, message: \n{command.stderr}")

    def _execute_tool(self) -> None:
        """
        Run a GATK function. Overrules method in parent class
        :return: None
        """
        self._set_input()
        self._set_output()
        self._set_specific_parameters()
        self._build_command()
        try:
            self._execute_command()
        finally:
            shutil.rmtree(self._temp_dir)
        self._print_output()
        self._set_informs()
