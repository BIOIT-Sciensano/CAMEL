import re
from typing import Optional

from camel.app.tools.picard.picard import Picard


class MergeBamAlignment(Picard):
    """
    Class for Picard MergeBamAlignment function

    Required inputs:
    ----------------
    'BAM_UNMAPPED':     ToolIOFile object. Input uBAM file.

    Optional inputs:
    ----------------
    'BAM_ALIGNED'       ToolIOFile object. Input BAM file
    'FASTA_REF':        ToolIOFile object. FASTA file containing the reference genome.

    Output:
    -------
    'BAM':              ToolIOFile object. Output BAM file containing recalibrated read data.
    """

    def __init__(self):
        """
        Initialize a picard tool
        :return: None
        """
        super().__init__('Picard MergeBamAlignment', '2.23.3')

        self._required_inputs = ['BAM_UNMAPPED']
        self._specific_parameters = ["attributes_to_remove_multi"]

    def _set_input(self) -> None:
        """
        Set required inputs specification. Adds the specific input files for this tool to the input set
        by the parent class
        :return: None
        """
        # Initialize input string in parent class
        super()._set_input()

        # Add specific MergeBamAlignment files
        self._input_string += f"UNMAPPED={self._tool_inputs['BAM_UNMAPPED'][0].path} "

        if 'BAM_ALIGNED' in self._tool_inputs:
            self._input_string +=  f"ALIGNED={self._tool_inputs['BAM_ALIGNED'][0].path} "

    def _build_command(self) -> None:
        """
        Build the command to run tool
        :return: None
        """
        build_options = self._build_options(excluded_parameters=self._specific_parameters, delimiter='=')

        # Attributes from the alignment record that should be removed when merging; multiple options possible
        if 'attributes_to_remove_multi' in self._parameters:
            build_options.append(self.__split_multi_options('attributes_to_remove_multi'))

        option_string = " ".join(build_options)

        self._command.command = " ".join([
            "java", self._java_options, "-jar $PICARD_JAR", self._tool_command, self._java_options_temp_dir,
            self._input_string, self._output_string, option_string
        ])

    def _set_informs(self, stderr: Optional[str] = None) -> None:
        """
        Analyse the result of picard run and update tool.informs
        :return: None
        """
        for line in (self._command.stderr if stderr is None else stderr).splitlines():
            m = re.search(r'Finished reading (\d+) total records from alignment SAM/BAM.', line)
            if m:
                self.informs['reads_total'] = m.group(1)
            m = re.search(r'Wrote (\d+) alignment records and (\d+) unmapped reads.', line)
            if m:
                if m.group(1):
                    self.informs['aligned_reads'] = m.group(1)
                if m.group(2):
                    self.informs['unmapped_reads'] = m.group(2)

    def __split_multi_options(self, option) -> str:
        """
        Multiple values allowed for certain parameters. These are passed in a comma separated string and need to be split
        """
        option_list = self._parameters[option].value.split(",")
        return "".join(f"{self._parameters[option].option}={s} " for s in option_list)
