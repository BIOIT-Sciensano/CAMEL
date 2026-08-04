from camel.app.tools.picard.picard import Picard


class AddOrReplaceReadGroups(Picard):
    """
    Class for Picard AddOrReplaceReadGroups function
    """

    def __init__(self):
        """
        Initialize a picard tool
        :return: None
        """
        super().__init__('Picard AddOrReplaceReadGroups', '2.23.3')

    def _set_input(self) -> None:
        """
        Set input for the picard function
        :return: None
        """
        super()._set_input()
        if 'SAMPLE_NAME' in self._tool_inputs:
            # if SAMPLE_NAME specified, it will replace the default values of parameters: RG_sample_name, RG_id in DB
            self._specific_parameters = ['RG_id', 'RG_sample_name']
            self._input_string += f" RGSM={self._tool_inputs['SAMPLE_NAME'][0].value} RGID={self._tool_inputs['SAMPLE_NAME'][0].value}"
