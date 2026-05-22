import re
from typing import Optional

from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.picard.picard import Picard


class MarkDuplicates(Picard):
    """
    Class for Picard MarkDuplicates function

    Required inputs:
    ----------------
    'BAM':      ToolIOFile object. Input BAM file.
    'FASTA_REF':ToolIOFile object. FASTA file containing the reference genome.

    Output:
    -------
    'BAM':              ToolIOFile object. Output BAM file duplicate marked data.
    'METRICS':          ToolIOFile object. Output txt file containing metrics.
    """

    def __init__(self):
        """
        Initialize a picard tool
        :return: None
        """
        super().__init__('Picard MarkDuplicates', '2.23.3')

    def _set_input(self) -> None:
        """
        Set the input specification. This method overrides the method in the parent class.
        :return: None
        """
        # input can be one or multiple BAM or SAM files
        self._input_string += ''.join(f"I={f.path} " for f in self._tool_inputs[self._main_input])

        # Remove so it's not set again by the parent class function
        del self._tool_inputs[self._main_input]

        # Run parent function for other inputs, e.g. FASTA_REF (optional)
        super()._set_input()

    def _set_informs(self, stderr: Optional[str] = None) -> None:
        """
        Analyse the result of picard run and update tool.informs
        Total reads: Total no. or reads in the file
        Pairs unmatched: Total no. of unmatched pairs
        Duplicates count: Total no. or records marked as duplicates
        Optical duplicate cluster count: Total no. of optical duplicate clusters
        :return: None
        """
        for line in (self._command.stderr if stderr is None else stderr).splitlines():
            m = re.search(r'Read (\d+) records. (\d+) pairs never matched', line)
            if m:
                self.informs['reads_total'] = m.group(1)
                self.informs['pairs_unmatched'] = m.group(2)

            m = re.search(r'Marking (\d+) records as duplicates', line)
            if m:
                self.informs['duplicates_count'] = m.group(1)

            m = re.search(r'Found (\d+) optical duplicate clusters', line)
            if m:
                self.informs['optical_duplicate_clusters_count'] = m.group(1)

    def _set_output(self) -> None:
        """
        Set the output specification. Overrides parent class method.
        :return: None
        """
        self._tool_outputs['BAM'] = [ToolIOFile(self.folder / self._parameters['output'].value)]
        self._tool_outputs['METRICS'] = [
            ToolIOFile(self.folder / self._parameters['metrics_output'].value)]
