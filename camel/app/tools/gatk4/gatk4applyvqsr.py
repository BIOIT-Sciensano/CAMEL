from camel.app.error import InvalidToolInputError
from camel.app.tools.gatk4.gatk4 import GATK4


class GATK4ApplyVQSR(GATK4):
    """
    =============================
    GATK ApplyVQSR 4.1.9.0
    =============================
    Apply a score cutoff to filter variants based on a recalibration table.

    Inputs:
    -------
    'TXT_RecalibrationTable':   ToolIOFile object. The input recal file used by ApplyRecalibration
    'VCF':                      ToolIOFile object. One or more VCF files containing variants
    'TXT_tranches' (optional):  ToolIOFile object. The tranches file

    Output:
    -------
    'VCF'       The output filtered and recalibrated VCF file in which each variant is annotated with its VQSLOD value

    """

    def __init__(self) -> None:
        """
        Initialize GATKApplyVQSR tool.
        :return: None
        """
        super().__init__('gatk4 ApplyVQSR', '4.1.9.0')
        self._required_inputs = ['VCF', 'TXT_RecalibrationTable']
        self._output_type = 'VCF'

    def _set_input(self) -> None:
        """
        Set GATKApplyVQSR input. Adds recalibration table to parent class input
        :return: None
        """
        super()._set_input()
        if 'TXT_RecalibrationTable' in self._tool_inputs:
            self._input_string += f"--recal-file {self._tool_inputs['TXT_RecalibrationTable'][0].path} "

        if 'filter_level' in self._parameters:
            if 'TXT_tranches' in self._tool_inputs:
                self._input_string += f"--tranches-file {self._tool_inputs['TXT_tranches'][0].path} "
            else:
                raise InvalidToolInputError("Filter level parameter requires a tranches file.")
