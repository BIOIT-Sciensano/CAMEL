from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.picard.picard import Picard


class IntervalListTools(Picard):
    """
    ===============================
    Picard IntervalListTools 2.23.3
    ===============================
    This tool offers multiple interval list file manipulation capabilities, including: sorting, merging, subtracting,
    padding, and other set-theoretic operations.

    Required inputs:
    ----------------
    ['TXT_intervals'|'VCF'] ToolIOFile object. One or more interval lists. Supported formats are interval_list and VCF.
                            IntervalList should be denoted with the extension .interval_list, while a VCF must have one
                            of .vcf, .vcf.gz, .bcf

    Output:
    -------
    'TXT_intervalListFiles' ToolIOFile object. One or more interval lists

    """
    def __init__(self):
        """
        Initialize a picard tool
        :return: None
        """
        super().__init__('Picard IntervalListTools', '2.23.3')
        self._required_inputs = ['TXT_intervals', 'VCF']
        self._output_type = 'TXT_intervalLists'

    def _check_input(self):
        """
        Check input for a tool and prepare command line parameters for input
        :return: None
        """
        #  One or more input files possible - can be of type 'interval_list' or 'VCF'
        if ('TXT_intervals' in self._tool_inputs) and ('VCF' in self._tool_inputs):
            raise InvalidToolInputError()
        elif 'VCF' in self._tool_inputs:
            self._required_inputs.remove('TXT_intervals')
        elif 'TXT_intervals' in self._tool_inputs:
            self._required_inputs.remove('VCF')

        super()._check_input()

    def _set_input(self) -> None:
        """
        Set the input specification
        Overrides method in parent class.
        :return: None
        """
        if 'TXT_intervals' in self._tool_inputs:
            for interval_file in self._tool_inputs['TXT_intervals']:
                self._input_string += f"INPUT={interval_file.path} "

        if 'VCF' in self._tool_inputs:
            for vcf_file in self._tool_inputs['VCF']:
                self._input_string += f"INPUT={vcf_file.path} "

    def _set_output(self) -> None:
        """
        Set the output specification. Overrides method in parent class.
        - Scatter count < 1 or not specified: ToolIOFile object containing one interval list output file
        - Scatter count > 1: A list of ToolIOFile objects, n=scatter_count; automated naming of interval list
          files and folder cfr. structure described in code
        :return: None
        """
        scatter_count = int(self._parameters['scatter_count'].value) if 'scatter_count' in self._parameters else 1

        if scatter_count == 1:
            self._tool_outputs['TXT_intervalLists'] = [ToolIOFile(self._folder / self._parameters['output'].value)]
        else:
            intervallists_list = []
            for n in range(1, scatter_count + 1):
                intervallists_list.append(ToolIOFile(self._folder / self._parameters['output'].value / f'temp_{n:04d}_of_{scatter_count:d}' / "scattered.interval_list"))
            self._tool_outputs['TXT_intervalLists'] = intervallists_list
