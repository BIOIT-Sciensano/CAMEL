from camel.app.tools.gatk4.gatk4 import GATK4


class GATK4GatherBQSRReports(GATK4):
    """
    =============================
    GATK GatherBQSRReports 4.1.9.0
    =============================
    Gathers scattered BQSR recalibration reports into a single file

    Required inputs:
    ----------------
    'TXT_intervals':            ToolIOFile object. List of scattered BQSR report files

    Output:
    -------
    'TXT_RecalibrationTable':   ToolIOFile object. Text file containing the gathered recalibration data.
    """

    def __init__(self) -> None:
        """
        Initialize GATK4GatherBQSRReports tool.
        :return: None
        """
        super().__init__('gatk4 GatherBQSRReports', '4.1.9.0')
        self._required_inputs = ['TXT_intervals']
        self._output_type = 'TXT_RecalibrationTable'

    def _set_input(self) -> None:
        """
        Set the input specification in the input_string
        Overrides method in parent class.
        :return: None
        """
        # set input reports
        self._input_string += ''.join(f" --input {f.path}" for f in self._tool_inputs['TXT_intervals'])
