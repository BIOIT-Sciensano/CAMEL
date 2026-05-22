import os
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class Itsx(Tool):
    """
    ITSx is an open source software utility to extract the highly variable ITS1 and ITS2 subregions from ITS sequences,
    which is commonly used as a molecular barcode for e.g. fungi. As the inclusion of parts of the neighbouring, very
    conserved, ribosomal genes (SSU, 5S and LSU rRNA sequences) in the sequence identification process can lead to
    severely misleading results, ITSx identifies and extracts only the ITS regions themselves.
    """

    def __init__(self) -> None:
        """
        Initialize tool
        :return: None
        """
        super().__init__('itsx', '1.1.3')

    def _execute_tool(self) -> None:
        """
        Runs ITSx
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid:
        - FASTA key is required
        - Only one input file allowed
        - No other input keys are allowed
        :return: None
        """
        super()._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError(f'Not enough valid input files given for ITSx (FASTA required): {self._tool_inputs!r}')
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidToolInputError(f'Invalid number (max = 1) of FASTA files given for ITSx: {self._tool_inputs!r}')
        if len(self._tool_inputs.keys()) > 2:
            raise InvalidToolInputError(f'Too many input keys given for ITSx (only FASTA allowed): {self._tool_inputs!r}')

    def __get_basename(self) -> str:
        """
        Returns the prefix that will be used in the output.
        :return: String with the prefix used in the output
        """
        infile = self._tool_inputs['FASTA'][0].basename
        return f'{os.path.join(self._folder, os.path.splitext(infile)[0])}_ITSx'

    def __set_output(self) -> None:
        """
        Sets the name of the output files
        :return: None
        """
        basename = self.__get_basename()
        self._tool_outputs['TEXT_Summary'] = [ToolIOFile(Path(basename + '.summary.txt'))]
        self._tool_outputs['TEXT_NoDetection'] = [ToolIOFile(Path(basename + '_no_detections.fasta'))]
        self._tool_outputs['TEXT_Problematic'] = [ToolIOFile(Path(basename + '.problematic.txt'))]
        self._tool_outputs['TSV_Positions'] = [ToolIOFile(Path(basename + '.positions.txt'))]
        self._tool_outputs['GRAPH'] = [ToolIOFile(Path(basename + '.graph'))]
        self._tool_outputs['FASTA_Full'] = [ToolIOFile(Path(basename + '.full.fasta'))]
        self._tool_outputs['FASTA_ITS1'] = [ToolIOFile(Path(basename + '.ITS1.fasta'))]
        self._tool_outputs['FASTA_ITS2'] = [ToolIOFile(Path(basename + '.ITS2.fasta'))]
        self._tool_outputs['FASTA_NoDetection'] = [ToolIOFile(Path(basename + '_no_detections.fasta'))]
        if os.path.isfile(basename + '.chimeric.fasta'):
            self._tool_outputs['FASTA_Chimeric'] = [ToolIOFile(Path(basename + '.chimeric.fasta'))]

    def __build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        return ' '.join([
            f"-i {self._tool_inputs['FASTA'][0].path}",
            f'-o {self.__get_basename()}'
        ])

    def __build_command(self) -> None:
        """
        Concatenates required parameters and options to build the command to run 'make.contigs'
        :return: None
        """
        input_string = self.__build_input_string()
        options_string = ' '.join(self._build_options())
        self._command.command = ' '.join([self._tool_command, input_string, options_string])
