import json

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool

DESCRIPTION_BY_RES_TYPE = {
    'NONE': 'No resistance',
    'MONO': 'Mono-resistance',
    'MDR': 'Multi-drug resistance (MDR)',
    'PRE_XDR': 'Pre-extensive drug resistance (pre-XDR)',
    'XDR': 'Extensive drug resistance (XDR)',
    'OTHER': 'Other'
}


class AMRTypeDetermination(Tool):
    """
    Determines the type of AMR based on the predicted resistances.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Mycobacterium: AMR type determination', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'JSON' not in self._tool_inputs:
            raise InvalidToolInputError('Predicted resistance input is required (JSON)')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Parse input data
        with open(self._tool_inputs['JSON'][0].path) as handle:
            ab_data_all = json.load(handle)

        # Check resistances
        is_first_line_resistant = all(
            ab_data['phenotype'].startswith('R') for ab_data in ab_data_all if ab_data['category'] == 'First line')
        is_second_line_resistant_a = any(
            ab_data['phenotype'].startswith('R') for ab_data in ab_data_all if ab_data['category'] ==
            'Second line (group A)')
        is_second_line_resistant_b = any(
            ab_data['phenotype'].startswith('R') for ab_data in ab_data_all if ab_data['category'] ==
            'Second line (group B)')
        nb_resistances = sum(ab_data['phenotype'].startswith('R') for ab_data in ab_data_all)

        # Determine resistance type
        if is_first_line_resistant and is_second_line_resistant_a and is_second_line_resistant_b:
            res_type = 'XDR'
        elif is_first_line_resistant and (is_second_line_resistant_a or is_second_line_resistant_b):
            res_type = 'PRE_XDR'
        elif is_first_line_resistant:
            res_type = 'MDR'
        elif nb_resistances == 1:
            res_type = 'MONO'
        elif nb_resistances == 0:
            res_type = 'NONE'
        else:
            res_type = 'OTHER'

        # Create output file
        output_file = self.folder / 'resistance_type.json'
        with output_file.open('w') as handle:
            json.dump({
                'first_line_resistant': is_first_line_resistant,
                'second_line_group_a_resistant': is_second_line_resistant_a,
                'second_line_group_b_resistant': is_second_line_resistant_b,
                'resistance_type': res_type,
                'resistance_type_full': DESCRIPTION_BY_RES_TYPE[res_type]
            }, handle, indent=2)
        self._tool_outputs['JSON'] = [ToolIOFile(output_file)]
