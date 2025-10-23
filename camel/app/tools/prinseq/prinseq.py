import os

from camel.app.core.command import Command
from camel.app.core.errors import ToolExecutionError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool


class Prinseq(Tool):
    """
    PRINSEQ can be used to filter, reformat, or trim your genomic and metagenomic sequence data.
    """

    def __init__(self):
        """
        Initialize tool
                :return: None
        """
        super().__init__('prinseq', '0.20.4')

    def _execute_tool(self):
        """
        Runs Prinseq
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTQ_PE, FASTQ_SE, FASTA_PE or FASTA_SE key is required
        - Only one input file allowed for FASTQ_SE or FASTA_SE
        - Only two input files allowed for FASTQ_PE or FASTA_SE
        - No other input keys are allowed
        :return: None
        """
        super(Prinseq, self)._check_input()
        if len(self._tool_inputs.keys()) != 1:
            raise ValueError('Too many input keys given voor PRINSEQ: {!r}'.format(self._tool_inputs))
        if list(self._tool_inputs.keys())[0] not in ['FASTQ_PE', 'FASTA_PE', 'FASTQ_SE', 'FASTA_SE']:
            raise ValueError('Invalid input key given for PRINSEQ: {!r}'.format(self._tool_inputs))
        key, value = list(self._tool_inputs.items())[0]
        if key.endswith('PE') and len(self._tool_inputs[key]) != 2:
            raise ValueError('Invalid number (!= 2) of files given for PE key'
                             'for PRINSEQ: {!r}'.format(self._tool_inputs))
        elif key.endswith('SE') and len(self._tool_inputs[key]) != 1:
            raise ValueError('Invalid number (!= 1) of files given for SE key'
                             'for PRINSEQ: {!r}'.format(self._tool_inputs))

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command to run
        :return: None
        """
        input_string = self.__build_input_string()
        options_string = ' '.join(self._build_options())
        self._command.command = ' '.join([self._tool_command, input_string, options_string])

    def __build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        # Only one key is allowed
        input_key = list(self._tool_inputs.keys())[0]
        basename = self.__get_basename()
        filetype = 'fastq' if input_key.startswith('FASTQ') else 'fasta'
        input_string = '-{} {}'.format(filetype, self._tool_inputs[input_key][0].path)
        if input_key.endswith('PE'):
            input_string += ' -{}2 {}'.format(filetype, self._tool_inputs[input_key][1].path)
        input_string += ' -out_good {}'.format(basename + '.prinseq')
        input_string += ' -out_bad null'
        return input_string

    def __get_basename(self):
        """
        Returns the prefix that will be used in the output.
        :return: String with the prefix used in the output
        """
        infile = os.path.basename(list(self._tool_inputs.values())[0][0].path)
        return os.path.join(self._folder, infile[:infile.rfind('.')])

    def _check_command_output(self, command: Command):
        """
        Checks if the command was executed successfully.
        :param command: Command to check
        :return: None
        """
        for line in command.stderr.splitlines():
            if 'ERROR' in line:
                raise ToolExecutionError(self.name, f"Command execution failed (stderr: {command.stderr}).")
        if self._command.exit_code != 0:
            raise ToolExecutionError(self.name, f"Command execution failed (Exit code: {command.exit_code})")

    def __set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        basename = self.__get_basename()
        # Only one key is allowed
        input_key = list(self._tool_inputs.keys())[0]
        format_key = input_key[:6]
        filetype = 'fastq' if input_key.startswith('FASTQ') else 'fasta'
        if input_key.endswith('PE'):
            self._tool_outputs[format_key + 'PE'] = [ToolIOFile(basename + '.prinseq_1.' + filetype),
                                                     ToolIOFile(basename + '.prinseq_2.' + filetype)]
            self._tool_outputs[format_key + 'SE'] = [ToolIOFile(basename + '.prinseq_1_singletons.' + filetype),
                                                     ToolIOFile(basename + '.prinseq_2_singletons.' + filetype)]
        else:
            self._tool_outputs[format_key + 'SE'] = [ToolIOFile(basename + '.prinseq.' + filetype)]
