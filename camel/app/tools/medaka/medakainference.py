from camel.app.tools.medaka.medaka import Medaka


class MedakaInference(Medaka):
    """
    Class for Medaka inference function.
    Runs the medaka inference algorithm and outputs a HDF file in preparation for building a consensus sequence.
    """

    def __init__(self) -> None:
        """
        Initializes Medaka inference.
        :return: None
        """
        super().__init__('medaka inference')

        self._required_inputs = ['BAM']
        self._output_type = 'HDF'

    def _set_input(self) -> None:
        """
        Sets the input specifications and the input string.
        :return: None
        """
        super()._set_input()
        self._input_string = str(self._tool_inputs['BAM'][0].path)
