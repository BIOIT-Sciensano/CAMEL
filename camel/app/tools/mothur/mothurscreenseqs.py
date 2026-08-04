from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.mothur.mothur import Mothur


class MothurScreenSeqs(Mothur):
    """
    The screen.seqs command enables you to keep sequences that fulfill certain user defined criteria. Furthermore, it
    enables you to cull those sequences not meeting the criteria from a names, group, contigsreport, alignreport and
    summary file.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_screen_seqs')
        self._required_input = ['FASTA']
        self._optional_input = [
            'TSV_Summary',
            'TSV_Groups',
            'TSV_AlignReport',
            'TSV_ContigReport',
            'TSV_Names',
            'TSV_Counts',
            'TSV_Qfile',
            'TSV_Taxonomy',
        ]

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and input/output directories
        Example: ffastq=FileR1.fastq, rfastq=FileR2.fastq, inputdir=/test/data/input/,
        outputdir=/test/data/outputdir
        :return: String with the input parameters
        """
        items = []
        # As there can be only one file per key, this first file of the list is added
        for key, input_files in self._tool_inputs.items():
            if key == 'FASTA':
                items.append(f'fasta={input_files[0]}')
            elif key == 'TSV_Groups':
                items.append(f'group={input_files[0]}')
            elif key == 'TSV_Summary':
                items.append(f'summary={input_files[0]}')
            elif key == 'TSV_Names':
                items.append(f'name={input_files[0]}')
            elif key == 'TSV_AlignReport':
                items.append(f'alignreport={input_files[0]}')
            elif key == 'TSV_ContigReport':
                items.append(f'contigreport={input_files[0]}')
            elif key == 'TSV_Taxonomy':
                items.append(f'taxonomy={input_files[0]}')
            elif key == 'TSV_Counts':
                items.append(f'count={input_files[0]}')
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        # Screen.seqs re-appends the original extension (e.g. .fasta or .align) after first adding .good
        fasta_extension = self._get_basename().suffix
        self._tool_outputs['FASTA'] = [
            ToolIOFile(self._get_basename().with_suffix(f'.good{fasta_extension}'))
        ]
        if self._get_basename().with_suffix('.bad.accnos').exists():
            self._tool_outputs['TSV_Bad'] = [
                ToolIOFile(self._get_basename().with_suffix('.bad.accnos'))
            ]
        # Depending on the input keys more files will be created
        if 'TSV_Groups' in self._tool_inputs:
            self._tool_outputs['TSV_Groups'] = [
                ToolIOFile(self._get_basename('TSV_Groups').with_suffix('.good.groups'))
            ]
        if 'TSV_Summary' in self._tool_inputs:
            self._tool_outputs['TSV_Summary'] = [
                ToolIOFile(
                    self._get_basename('TSV_Summary').with_suffix('.good.summary')
                )
            ]
        if 'TSV_Names' in self._tool_inputs:
            self._tool_outputs['TSV_Names'] = [
                ToolIOFile(self._get_basename('TSV_Names').with_suffix('.good.names'))
            ]
        if 'TSV_AlignReport' in self._tool_inputs:
            self._tool_outputs['TSV_AlignReport'] = [
                ToolIOFile(
                    self._get_basename('TSV_AlignReport', {'.report'}).with_suffix(
                        '.good.align.report'
                    )
                )
            ]
        if 'TSV_ContigReport' in self._tool_inputs:
            self._tool_outputs['TSV_ContigReport'] = [
                ToolIOFile(
                    self._get_basename('TSV_ContigReport', {'.report'}).with_suffix(
                        '.good.contigs.report'
                    )
                )
            ]
        if 'TSV_Counts' in self._tool_inputs:
            self._tool_outputs['TSV_Counts'] = [
                ToolIOFile(
                    self._get_basename('TSV_Counts').with_suffix('.good.count_table')
                )
            ]
