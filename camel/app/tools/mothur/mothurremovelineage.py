from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.mothur.mothur import Mothur


class MothurRemoveLineage(Mothur):
    """
    The remove.lineage command reads a taxonomy file and a taxon and generates a new file that contains only the
    sequences not containing that taxon.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_remove_lineage')
        self._required_input = ['TSV_Taxonomy']
        self._optional_input = [
            'FASTA',
            'TSV_Names',
            'TSV_Counts',
            'TSV_Groups',
            'TSV_AlignReport',
            'TSV_List',
        ]

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = [f"taxonomy={self._tool_inputs['TSV_Taxonomy'][0]}"]
        input_parameters = {
            'FASTA': 'fasta=',
            'TSV_Names': 'name=',
            'TSV_Counts': 'count=',
            'TSV_Groups': 'group=',
            'TSV_AlignReport': 'alignreport=',
            'TSV_List': 'list=',
        }
        for key, input_files in self._tool_inputs.items():
            # Based on the key the correct option flag is added to the input string
            if key != 'TSV_Taxonomy':
                items.append(f'{input_parameters[key]}{input_files[0].path}')
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        # The extension of the output files is not always built in the same way. A dictionary
        # holds the suffixes that can be used when getting the base name
        output_extensions = {
            'FASTA': [None, '.pick.fasta'],
            'TSV_Names': [None, '.pick.names'],
            'TSV_Counts': [None, '.pick.count_table'],
            'TSV_Groups': [None, '.pick.groups'],
            'TSV_AlignReport': [{'.report'}, '.pick.align.report'],
            'TSV_List': [None, '.pick.list'],
            'TSV_Taxonomy': [None, '.pick.taxonomy'],
        }
        for key, input_files in self._tool_inputs.items():
            basename = self._get_basename(key, output_extensions[key][0])
            self._tool_outputs[key] = [
                ToolIOFile(basename.with_suffix(output_extensions[key][1]))
            ]
        self._tool_outputs['TSV_Accnos'] = [
            ToolIOFile(self._get_basename('TSV_Taxonomy').with_suffix('.accnos'))
        ]
