from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothur import Mothur


class MothurScreenSeqs(Mothur):
    """
    The screen.seqs command enables you to keep sequences that fulfill certain user defined criteria. Furthermore, it
    enables you to cull those sequences not meeting the criteria from a names, group, contigsreport, alignreport and
    summary file.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_screen_seqs', version=None)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA key is required
        - Only 'FASTA', 'TSV_Summary', 'TSV_Groups', 'TSV_AlignReport', 'TSV_ContigReport',
          'TSV_Names', 'TSV_Counts', 'TSV_Qfile', and 'TSV_Taxonomy' are allowed
        - Only one input file per key is allowed
        :return: None
        """
        super()._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError(f'No input file given for Mothur screen.seqs: {self._tool_inputs!r}')
        for key, input_files in self._tool_inputs.items():
            if key not in ['FASTA', 'TSV_Summary', 'TSV_Groups', 'TSV_AlignReport', 'TSV_ContigReport',
                           'TSV_Names', 'TSV_Counts', 'TSV_Qfile', 'TSV_Taxonomy']:
                raise InvalidToolInputError(f'Invalid input key given for Mothur screen.seqs: {self._tool_inputs!r}')
            if len(input_files) != 1:
                raise InvalidToolInputError(f'Invalid number (max = 1) of files given for Mothur \
                                                     screen.seqs: {self._tool_inputs!r}')

    def _build_input_string(self):
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
            elif key == 'TSV_ContigsReport':
                items.append(f'contigreport={input_files[0]}')
            elif key == 'TSV_Taxonomy':
                items.append(f'taxonomy={input_files[0]}')
            elif key == 'TSV_Counts':
                items.append(f'count={input_files[0]}')
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        # Screen.seqs re-appends the original extension (e.g. .fasta or .align) after first adding .good
        path = self._tool_inputs['FASTA'][0].path
        fasta_extension = path[self._tool_inputs['FASTA'][0].path.rfind('.'):]
        self._tool_outputs['FASTA'] = [ToolIOFile(super()._get_basename('FASTA') + '.good' + fasta_extension)]
        self._tool_outputs['TSV_Bad'] = [ToolIOFile(super()._get_basename('FASTA') + '.bad.accnos')]
        # Depending on the input keys more files will be created
        if 'TSV_Groups' in self._tool_inputs:
            self._tool_outputs['TSV_Groups'] = [ToolIOFile(super()._get_basename('TSV_Groups', '.groups') +
                                                '.good.groups')]
        if 'TSV_Summary' in self._tool_inputs:
            self._tool_outputs['TSV_Summary'] = [ToolIOFile(super()._get_basename('TSV_Summary', '.summary') +
                                                 '.good.summary')]
        if 'TSV_Names' in self._tool_inputs:
            self._tool_outputs['TSV_Names'] = [ToolIOFile(super()._get_basename('TSV_Names', '.names') +
                                               '.good.names')]
        if 'TSV_AlignReport' in self._tool_inputs:
            self._tool_outputs['TSV_AlignReport'] = [ToolIOFile(super()._get_basename('TSV_AlignReport', '.align.report') +
                                                     '.good.align.report')]
        if 'TSV_ContigsReport' in self._tool_inputs:
            self._tool_outputs['TSV_ContigsReport'] = [ToolIOFile(super()._get_basename('TSV_ContigReport', '.contig.report') +
                                                       '.good.contigs.report')]
        if 'TSV_Counts' in self._tool_inputs:
            self._tool_outputs['TSV_Counts'] = [ToolIOFile(super()._get_basename('TSV_Counts', '.count_table') +
                                                '.good.count_table')]
