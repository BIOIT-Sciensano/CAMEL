from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.mothur.mothur import Mothur


class MothurClassifyOtu(Mothur):
    """
    The classify.otu command is used to get a consensus taxonomy for an otu.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_classify_otu')
        self._required_input = ['TSV_List', 'TSV_Taxonomy']
        self._optional_input = [
            'TSV_Groups',
            'TSV_Counts',
            'TSV_Names',
            'TSV_RefTaxonomy',
        ]

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = [
            f"list={self._tool_inputs['TSV_List'][0]}",
            f"taxonomy={self._tool_inputs['TSV_Taxonomy'][0]}",
        ]
        if 'TSV_Counts' in self._tool_inputs:
            items.append(f"count={self._tool_inputs['TSV_Counts'][0]}")
        if 'TSV_Groups' in self._tool_inputs:
            items.append(f"group={self._tool_inputs['TSV_Groups'][0]}")
        if 'TSV_Names' in self._tool_inputs:
            items.append(f"name={self._tool_inputs['TSV_Names'][0]}")
        if 'TSV_RefTaxonomy' in self._tool_inputs:
            items.append(f"reftaxonomy={self._tool_inputs['TSV_RefTaxonomy'][0]}")
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them.
        :return: None
        """
        labels = self._get_labels()
        self._tool_outputs.update({'TSV_Taxonomy': [], 'TSV_Summary': []})
        basename = self._get_basename('TSV_List')
        # Each label creates a seperate output
        for label in labels:
            self._tool_outputs['TSV_Taxonomy'] += [
                ToolIOFile(basename.with_suffix(f'.{label}.cons.taxonomy'))
            ]
            self._tool_outputs['TSV_Summary'] += [
                ToolIOFile(basename.with_suffix(f'.{label}.cons.tax.summary'))
            ]
