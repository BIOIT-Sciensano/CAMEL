from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class Tabix(Tool):
    """
    Tabix indexes a TAB-delimited genome position file in.tab.bgz and creates an index file (in.tab.bgz.tbi or
    in.tab.bgz.csi) when region is absent from the command-line. The input data file must be position sorted and
    compressed by bgzip which has a gzip(1) like interface.

    After indexing, tabix is able to quickly retrieve data lines overlapping regions specified in the format
    "chr:beginPos-endPos". (Coordinates specified in this region format are 1-based and inclusive.)
    """

    def __init__(self):
        """
        Initializes this tool.
        """
        super().__init__('tabix', '1.9')

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        output_filename = self._folder / self.get_param_value('output_filename')
        self._tool_outputs['TSV'] = [ToolIOFile(output_filename)]

    def _check_input(self):
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'TAB' not in self._tool_inputs:
            raise InvalidToolInputError("TAB index is required")
        if not any(key in self._tool_inputs for key in ['BED', 'VAL_regions']):
            raise InvalidToolInputError("Either BED or regions input is required")
        super()._check_input()

    def __build_command(self):
        """
        Builds the command line call.
        :return: None
        """
        parts = [
            self._tool_command,
            self._tool_inputs['TAB'][0].path,
        ]
        if 'BED' in self._tool_inputs:
            parts.insert(1, '--preset bed')
            parts.append('--regions {}'.format(self._tool_inputs['BED'][0].path))
        if 'VAL_regions' in self._tool_inputs:
            parts.extend([x.value for x in self._tool_inputs['VAL_regions']])
        parts.append('> {}'.format(self._parameters['output_filename'].value))
        self._command.command = ' '.join(parts)

    def __parse_output(self):
        """
        Parses the command output.
        :return: None
        """
        output = []
        for line in self._command.stdout.splitlines():
            output.append(line.strip().split('\t'))
        return output
