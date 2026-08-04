from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidParameterError, InvalidToolInputError
from camel.app.tools.bedtools.bedtools import Bedtools


class BedtoolsGenomecov(Bedtools):
    """
    Bedtools genomecov computes histograms (default), per-base reports (-d) and BEDGRAPH (-bg) summaries of feature
    coverage (e.g., aligned sequences) for a given genome.
    """
    OUTPUT_FILE_BASENAME = 'bedtools_genomecov'

    # Mutually exclusive, only ONE can be specified.
    OUTPUT_FORMAT_OPTIONS = ['BedGraphWithZeroCoverage', 'BedGraph', 'DepthWithZeroCoverage', 'Depth']

    def __init__(self) -> None:
        """
        Initialize a samtools tool.
        :return: None
        """
        super().__init__('bedtools genomecov', '2.31.0')
        self._output_filename = None
        self._required_inputs = ['BAM']

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        self._check_required_inputs()
        if len(self._tool_inputs['BAM']) != 1:
            raise InvalidToolInputError("Exactly one BAM input file expected.")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        params_fmt_out = next(param for param in self._parameters if param in BedtoolsGenomecov.OUTPUT_FORMAT_OPTIONS)
        extension = 'bed' if params_fmt_out.startswith('Bed') else 'tsv'
        path_out = self.folder / f"{BedtoolsGenomecov.OUTPUT_FILE_BASENAME}.{extension}"
        self.__build_command(path_out)
        self._execute_command()
        self._tool_outputs[extension.upper()] = [ToolIOFile(path_out)]

    def __build_command(self, path_out: Path) -> None:
        """
        Builds the command line call.
        :param path_out: Path to output file
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            *self._build_options(),
            f" -ibam {self._tool_inputs['BAM'][0].path}",
            f' > {path_out}'
        ])

    def _check_parameters(self) -> None:
        """
        Check the parameters
        :return: None
        """
        super()._check_parameters()
        params_excluded = []
        params_excluded += BedtoolsGenomecov.OUTPUT_FORMAT_OPTIONS
        output_opt_found = False

        for param in self._parameters.keys():
            if param in params_excluded:
                if not output_opt_found:
                    params_excluded.remove(param)
                    output_opt_found = True
                else:
                    raise InvalidParameterError("Only one output option should be specified.")
