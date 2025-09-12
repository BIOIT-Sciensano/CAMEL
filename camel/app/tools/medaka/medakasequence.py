from camel.app.tools.medaka.medaka import Medaka


class MedakaSequence(Medaka):
    """
    Class for Medaka sequence function.

    Medaka sequence reads the output of Medaka inference and outputs a consensus fasta sequence.
    """

    def __init__(self) -> None:
        """
        Initializes Medaka sequence.
        :return: None
        """
        super().__init__('medaka sequence', '2.0.0')

        self._required_inputs = ['HDF', 'FASTA']
        self._output_type = 'FASTA'

    def _set_input(self) -> None:
        """
        Sets the input specification and the input string.
        :return: None
        """
        super()._set_input()
        hdf_file = self._tool_inputs['HDF'][0].path
        fasta_file = self._tool_inputs['FASTA'][0].path
        self._input_string = f'{hdf_file} {fasta_file}'
