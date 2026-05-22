from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.mothur.mothur import Mothur


class MothurSummarySeqs(Mothur):
    """
    The summary.seqs command will summarize the quality of sequences in
    an unaligned or aligned fasta-formatted sequence file.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_summary_seqs')
        self._required_input = ['FASTA']
        self._optional_input = ['TSV_Counts']

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and input/output directories
        Example: fasta=File1.trim.contig.fasta, inputdir=/test/data/input/, outputdir=/test/data/outputdir
        :return: String with the input parameters
        """
        items = [f"fasta={self._tool_inputs['FASTA'][0]}"]
        if 'TSV_Counts' in self._tool_inputs:
            items.append(f"count={self._tool_inputs['TSV_Counts'][0]}")
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the output file object with them
        :return: None
        """
        output_base = self._folder / self._tool_inputs['FASTA'][0].basename
        self._tool_outputs['TSV_Summary'] = [ToolIOFile(output_base.with_suffix('.summary'))]
        self._tool_outputs['TSV_Stats'] = [ToolIOFile(output_base.with_suffix('.stats'))]

    def _execute_tool(self) -> None:
        """
        Runs Mothur summary.seqs
        :return: None
        """
        self._create_symlinks(self._temp_dir)
        self._build_command()
        self._execute_command()
        self.__write_stats_to_file()
        self._set_output()
        self.__set_informs()

    def __write_stats_to_file(self) -> None:
        """
        Writes the statistics that were output to stdout to a file
        :return: None
        """
        stats_file = self._get_basename().with_suffix('.stats')
        with open(stats_file, 'w', encoding='utf-8') as outhandle:
            write_to_output = False
            for line in self._command.stdout.splitlines():
                # The first line of the stats starts with two tabs (i.e. \t\tStart\tEnd...)
                if line.startswith('\t'):
                    write_to_output = True
                if write_to_output is True:
                    outhandle.write(line + '\n')
                # The last line of the stats starts with a '#' (i.e. # of Seqs...)
                if line.strip() == '' and write_to_output is True:
                    break

    def __set_informs(self) -> None:
        """
        Adds the summary statistics to the informs.
        :return: None
        """
        columns = ['Description', 'Start', 'End', 'NBases', 'Ambigs', 'Polymer', 'NumSeqs']
        with open(self._tool_outputs['TSV_Stats'][0].path, encoding='utf-8') as statsfile:
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
