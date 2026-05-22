from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothur import Mothur


class MothurChimeraUchime(Mothur):
    """
    The chimera.uchime command reads a fasta file and reference file and outputs potentially chimeric sequences.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_chimera_uchime')
        self._required_input = ['FASTA']
        self._optional_input = ['FASTA_Ref', 'TSV_Groups']

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid:
        - FASTA key is required
        - Either TSV_Names or TSV_Counts is required
        - FASTA_Ref and TSV_Groups are allowed as additional input
        - Only one input file per key is allowed
        - The use of TSV_Names is not yet implemented (lack of documentation)
        :return: None
        """
        if 'TSV_Counts' in self._tool_inputs:
            self._required_input.append('TSV_Counts')
        elif 'TSV_Names' in self._tool_inputs:
            self._required_input.append('TSV_Names')
        else:
            raise InvalidToolInputError('Either TSV_Counts or TSV_Names is required')
        super()._check_input()

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = [f"fasta={self._tool_inputs['FASTA'][0]}"]
        if 'TSV_Counts' in self._tool_inputs:
            items.append(f"count={self._tool_inputs['TSV_Counts'][0]}")
        elif 'TSV_Names' in self._tool_inputs:
            items.append(f"name={self._tool_inputs['TSV_Names'][0]}")
        if 'TSV_Groups' in self._tool_inputs:
            items.append(f"group={self._tool_inputs['TSV_Groups'][0]}")
        if 'FASTA_Ref' in self._tool_inputs:
            items.append(f"reference={self._tool_inputs['FASTA_Ref'][0]}")
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them.
        When dereplicate=false (default), mothur internally runs remove.seqs and the cleaned
        FASTA/count outputs carry a .pick prefix. When dereplicate=true, mothur only modifies
        the count table and keeps the .denovo.uchime prefix for all outputs.
        :return: None
        """
        basename = self._get_basename()
        self._tool_outputs['TSV_Chimeras'] = [
            ToolIOFile(basename.with_suffix('.denovo.uchime.chimeras'))
        ]
        self._tool_outputs['TSV_Accnos'] = [
            ToolIOFile(basename.with_suffix('.denovo.uchime.accnos'))
        ]
        dereplicate = ('dereplicate' in self._parameters and self.get_param_value('dereplicate') == 'true')
        if dereplicate:
            self._tool_outputs['FASTA'] = [
                ToolIOFile(basename.with_suffix('.denovo.uchime.fasta'))
            ]
            if 'TSV_Counts' in self._tool_inputs:
                self._tool_outputs['TSV_Counts'] = [
                    ToolIOFile(basename.with_suffix('.denovo.uchime.count_table'))
                ]
        else:
            self._tool_outputs['FASTA'] = [
                ToolIOFile(basename.with_suffix('.pick.fasta'))
            ]
            if 'TSV_Counts' in self._tool_inputs:
                self._tool_outputs['TSV_Counts'] = [
                    ToolIOFile(basename.with_suffix('.pick.count_table'))
                ]
        if 'TSV_Names' in self._tool_inputs:
            raise RuntimeError(
                'The use of a names file is not yet implemented for chimera.uchime as the '
                'outputs are unknown!'
            )
