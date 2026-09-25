from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidParameterError
from camel.app.tools.gatk.gatk import GATK


class GATKAnalyzeCovariates(GATK):
    """
    Generates plots and data about the covariates used during base quality score recalibration (BQSR).

    Required inputs:
    - TXT_TABLE_BEFORE: Recalibration table of the original data
    - TXT_TABLE_AFTER: Recalibration table of the recalibrated data
    - FASTA_REF: Reference genome

    Optional inputs:
    - BQSR: Additional recalibration table

    Outputs (at least one is required):
    - PDF: Recalibration plots (if 'pdf_output' is set, requires the R packages ggplot2, gplots, reshape and gsalib)
    - CSV: Recalibration covariates data (if 'csv_output' is set)
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('gatk AnalyzeCovariates')

    def _check_parameters(self) -> None:
        """
        Checks if the parameters are valid.
        :return: None
        """
        super()._check_parameters()
        if self.get_param_value('pdf_output') is None and self.get_param_value('csv_output') is None:
            raise InvalidParameterError("At least one output is required ('pdf_output' and/or 'csv_output')")

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['TXT_TABLE_BEFORE', 'TXT_TABLE_AFTER', 'FASTA_REF'])
        super()._check_input()

    def _get_input_options(self) -> list[str]:
        """
        Returns the command line options for the input files.
        :return: Input options
        """
        options = [
            *super()._get_input_options(),
            f"-before {self._tool_inputs['TXT_TABLE_BEFORE'][0].path}",
            f"-after {self._tool_inputs['TXT_TABLE_AFTER'][0].path}"
        ]
        if 'BQSR' in self._tool_inputs:
            options.append(f"-BQSR {self._tool_inputs['BQSR'][0].path}")
        return options

    def _set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        if self.get_param_value('pdf_output') is not None:
            self._tool_outputs['PDF'] = [ToolIOFile(self.folder / self.get_param_value('pdf_output'))]
        if self.get_param_value('csv_output') is not None:
            self._tool_outputs['CSV'] = [ToolIOFile(self.folder / self.get_param_value('csv_output'))]
