from camelcore.app.io.tooliovalue import ToolIOValue

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class SerotypeDetectorEcoli(Tool):
    """
    This tool detects the E. coli serotype based on gene hits.
    """

    def __init__(self):
        """
        Initializes this tool.
        """
        super().__init__('Serotype Detector', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        for char in ('H', 'O'):
            if f'HITS_{char}' not in self._tool_inputs:
                raise InvalidToolInputError(f"{char}-type hits are required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        h_type = self.__get_h_type()
        o_type = self.__get_o_type()
        self._tool_outputs['VAL_H_type'] = [ToolIOValue(h_type)]
        self._tool_outputs['VAL_O_type'] = [ToolIOValue(o_type)]
        self._tool_outputs['VAL_serotype'] = [ToolIOValue(f'{o_type}:{h_type}')]

    def __get_h_type(self) -> str:
        """
        Returns the H type based on the gene detection hits.
        :return: H type.
        """
        detected_genes = [io.value for io in self._tool_inputs['HITS_H']]
        if len(detected_genes) == 0:
            return 'H-unknown'
        elif len(detected_genes) == 1:
            return detected_genes[0].get_metadata_value('Predicted serotype')
        else:
            hits_sorted = {
                'fliC': [h for h in detected_genes if 'fliC' in h.locus],
                'non-fliC': [h for h in detected_genes if 'fliC' not in h.locus]}
            if len(hits_sorted['non-fliC']) == 1:
                return hits_sorted['non-fliC'][0].get_metadata_value('Predicted serotype')
        return 'H-ambiguous'

    def __get_o_type(self) -> str:
        """
        Returns the O type.
        :return: O type.
        """
        detected_genes = [io.value for io in self._tool_inputs['HITS_O']]
        if len(detected_genes) == 0:
            return 'O-unknown'
        elif len(detected_genes) == 1:
            return detected_genes[0].get_metadata_value('Predicted serotype')
        else:
            hits_sorted = {'xy': [h for h in detected_genes if any([g in h.locus for g in ('wzx', 'wzy')])],
                           'mt': [h for h in detected_genes if any([g in h.locus for g in ('wzm', 'wzt')])]}
            if len(hits_sorted['xy']) > 0 and len(hits_sorted['mt']) > 0:
                return 'O-ambiguous'
            else:
                all_hits = hits_sorted['xy'] + hits_sorted['mt']
                o_type = all_hits[0].get_metadata_value('Predicted serotype')
                for hit in all_hits:
                    if hit.get_metadata_value('Predicted serotype') != o_type:
                        return "O-ambiguous"
                return o_type
