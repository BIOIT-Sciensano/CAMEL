from pathlib import Path

from camel.app.components.files.fileutils import FileUtils
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.bwa.bwa import BWA


class BWAIndex(BWA):
    """
    Creates a BWA index for a FASTA file.
    """

    MULTI_FASTA_GENOME_FILE = 'complete_genome.fasta'

    def __init__(self) -> None:
        """
        Initialize BWA index
            :return: None
        """
        super().__init__('bwa_index', '0.7.17')
        self._refgenome_fasta = None

    def _execute_tool(self) -> None:
        """
        Function to run BWA index
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __get_multi_fasta_genome_filename(self) -> Path:
        """
        Get the filename used for multi fasta file representing complete genome
        :return: name of the multi fasta file with complete path
        """
        return self.folder / BWAIndex.MULTI_FASTA_GENOME_FILE

    def _check_input(self) -> None:
        """
        Check FASTA_REF input and concatenate them if multiple fasta input files
        :return: None
        """
        super()._check_input()

        nb_of_inputs = len(self._tool_inputs['FASTA_REF'])
        if nb_of_inputs == 0:
            raise ValueError("Required reference genome (FASTA) input file is missing.")
        elif nb_of_inputs > 1:
            multifasta_file = self.__get_multi_fasta_genome_filename()
            FileUtils.concatenate_files(Path(multifasta_file), [f.path for f in self._tool_inputs['FASTA_REF']])
            self._refgenome_fasta = multifasta_file
        else:
            self._refgenome_fasta = self.folder / self._tool_inputs['FASTA_REF'][0].basename
            if self._refgenome_fasta != self._tool_inputs['FASTA_REF'][0].path and not self._refgenome_fasta.exists():
                self._refgenome_fasta.symlink_to(self._tool_inputs['FASTA_REF'][0].path)

    def __set_output(self) -> None:
        """
        Set output for BWA index
        :return: None
        """
        self._tool_outputs['INDEX_GENOME_PREFIX'] = [ToolIOValue(self._refgenome_fasta)]

    def __build_command(self) -> None:
        """
        Build the command to run BWA index
        :return: None
        """
        options = ' '.join(self._build_options())
        self._command.command = f'{self._tool_command} {options} {self._refgenome_fasta}'
