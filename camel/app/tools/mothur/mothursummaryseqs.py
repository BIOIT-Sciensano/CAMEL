from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurSummarySeqs(Mothur):
    """
    The summary.seqs command will summarize the quality of sequences in
    an unaligned or aligned fasta-formatted sequence file.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_summary_seqs', '1.39.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA key is required
        - Allowed keys are 'FASTA', 'TSV_Counts'
        - Only one input file allowed per key
        :return: None
        """
        super()._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('No input file given for Mothur summary.seqs: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.items():
            if key not in ['FASTA', 'TSV_Counts']:
                raise InvalidToolInputError('Invalid input key given for Mothur summary.seqs: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidToolInputError('Invalid number (max = 1) of files given for Mothur \
                                                     summary.seqs: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and input/output directories
        Example: fasta=File1.trim.contig.fasta, inputdir=/test/data/input/,
        outputdir=/test/data/outputdir
        :return: String with the input parameters
        """
        items = ['fasta={}'.format(self._tool_inputs['FASTA'][0])]
        if 'TSV_Counts' in self._tool_inputs:
            items.append('count={}'.format(self._tool_inputs['TSV_Counts'][0]))
        items.append('outputdir={}'.format(self._folder))
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the output file object with them
        :return: None
        """
        basename = super()._get_basename()
        self._tool_outputs['TSV_Summary'] = [ToolIOFile(basename + '.summary')]
        self._tool_outputs['TSV_Stats'] = [ToolIOFile(basename + '.stats')]

    def _execute_tool(self):
        """
        Runs Mothur summary.seqs
        :return: None
        """
        self._create_symlinks()
        self._build_command()
        self._execute_command()
        self.__write_stats_to_file()
        self._symlink_cleanup()
        self._set_output()
        self.__set_informs()

    def __write_stats_to_file(self):
        """
        Writes the statistics that were output to stdout to a file
        :return: None
        """
        output_file = open(self._get_basename() + '.stats', 'wt', encoding='utf-8')
        write_to_output = False
        for line in self._command.stdout.splitlines():
            # The first line of the stats starts with two tabs (i.e. \t\tStart\tEnd...)
            if line.startswith('\t'):
                write_to_output = True
            if write_to_output is True:
                output_file.write(line + '\n')
            # The last line of the stats starts with a '#' (i.e. # of Seqs...)
            if line.strip() == '' and write_to_output is True:
                break

    def __set_informs(self):
        """
        Adds the summary statistics to the informs.
        :return: None
        """
        columns = ['Description', 'Start', 'End', 'NBases', 'Ambigs', 'Polymer', 'NumSeqs']
        with open(self._get_basename() + '.stats', 'rt', encoding='utf-8') as statsfile:
            for line in statsfile:
                if not line.startswith('\t\t'):
                    line_informs = line.strip().split('\t')
                    category = line_informs[0].strip()
                    if category in {'Minimum:', '2.5%-tile:', '25%-tile:', 'Median:', '75%-tile:', '97.5%-tile:', 'Maximum:'}:
                        self._informs[category[:-1]] = {}
                        for i in range(1, len(line_informs)):
                            self._informs[category[:-1]][columns[i]] = int(line_informs[i])
                    elif category == 'Mean:':
                        self._informs[category[:-1]] = {}
                        for i in range(1, len(line_informs)):
                            self._informs[category[:-1]][columns[i]] = float(line_informs[i])
                    elif category.lower().startswith('# of unique'):
                        self._informs['unique'] = int(line_informs[1])
                    elif category.lower().startswith('# of seqs') or category.lower().startswith('total # of seqs'):
                        self._informs['total'] = int(line_informs[1])
