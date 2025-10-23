import os
from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool


class Carma(Tool):
    """
    # TODO: deprecate
    CARMA is a software pipeline for the characterisation of species composition and the genetic potential of microbial
    samples using short reads. In contrast to the traditional 16S-rRNA approach for taxonomical classification, CARMA
    uses reads that encode for known proteins. By assigning the taxonomic origins to each read, a profile is
    constructed which characterises the taxonomic composition of the corresponding community.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('carma', '20150505')
        self._input_key = None
        self._carma_dir = 'carma3/20150505'
        self._hmmer_dir = 'hmmer/3.1b2'
        self._cmd_param = ['hmmer', 'blast', 'classify_egt', 'classify_blast', 'classify_rdp', 'type', 'blast_egts',
                           'fasta_input', 'database', 'delay']
        self._config_loc_param = ['blast_exe', 'blast_nr', 'blast_nt', 'ncbi', 'pfam']
        self._config_param = ['blastn_evalue', 'blastp_evalue', 'blastx_evalue', 'carma_blastn_alignment_length',
                              'carma_blastn_biscore', 'carma_blastn_evalue', 'carma_blastp_alignment_length',
                              'carma_blastp_bitscore', 'carma_blastp_evalue', 'carma_blastx_alignment_length',
                              'carma_blastx_bitscore', 'carma_blastx_evalue', 'clusteremail', 'hmmscan_evalue',
                              'lca_top_percent', 'max_blast_description_length', 'max_number_of_chunks', 'message',
                              'ncbi_max', 'pairwise_blosum_minoverlap', 'pairwise_blosum_minscore', 'query_multiplier',
                              'score_gapextension', 'score_gapopen', 'score_match', 'score_mismatch',
                              'use_hard_threshold', 'blosum_file']

    def _execute_tool(self):
        """
        Runs Carma
        :return: None
        """
        self.__set_input_key()
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - Either FASTA, BLASTX or EGT allowed
        - Only one input file allowed
        - No other input keys are allowed
        :return: None
        """
        super()._check_input()
        if [x in self._tool_inputs for x in ('FASTA', 'BLASTX', 'EGT')].count(True) != 1:
            raise InvalidToolInputError(f'Invalid input keys given for CARMA: {self._tool_inputs!r}')
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidToolInputError(f'Too many input keys given for CARMA: {self._tool_inputs!r}')
        for value in self._tool_inputs.values():
            if len(value) > 1:
                raise InvalidToolInputError(
                    f'Invalid number (max = 1) of files per key given for CARMA: {self._tool_inputs!r}')

    def _check_command_output(self, command: Command):
        """
        Checks if the command was executed successfully.
        :param command: Command to check
        :return: None
        """
        if 'ERROR' in self._command.stderr:
            raise ToolExecutionError(self.name, f"Command execution failed (stderr: {command.stderr}).")
        if self._command.exit_code != 0:
            raise ToolExecutionError(self.name, f"Command execution failed (Exit code: {command.exit_code})")

    def __set_input_key(self):
        """
        Sets the input key that is provided as input
        :return: None
        """
        self._input_key = list(self._tool_inputs.keys())[0]

    def __create_config_file(self):
        """
        Creates the config file that is used by CARMA from the values available in the database.
        :return: Location of the config file
        """
        option_list = super()._build_options(self._cmd_param + self._config_loc_param, ' = ')
        option_list += self.__build_fixed_options()
        option_list += self.__build_ncbi_options()
        option_list += self.__build_blast_options()
        option_list += self.__build_pfam_options()
        option_list += self.__build_unused_options()
        with open(os.path.join(self._folder, 'carma.cfg')) as handle:
            for option in option_list:
                handle.write(f'{option}\n')
        return os.path.join(self._folder, 'carma.cfg')

    def __build_ncbi_options(self):
        """
        Builds the options that specify the location of the NCBI database
        :return: List of options
        """
        ncbi = self._parameters['ncbi']
        return [
            f'nodes_dmp = /data/taxonomy/ncbi/{ncbi.value}/nodes.dmp',
            f'merged_dmp = /data/taxonomy/ncbi/{ncbi.value}/merged.dmp',
            f'names_dmp = /data/taxonomy/ncbi/{ncbi.value}/names.dmp'
        ]

    def __build_blast_options(self):
        """
        Builds the options that specify the location of the BLAST executables
        :return: List of options
        """
        blast = self._parameters['blast_exe']
        blast_nt = self._parameters['blast_nt']
        blast_nr = self._parameters['blast_nr']
        return [
            f'blastall_script = /usr/local/bin/blast/{blast.value}/blastall',
            f'fastacmd_script = /usr/local/bin/blast/{blast.value}/fastacmd',
            f'formatdb_script = /usr/local/bin/blast/{blast.value}/formatdb',
            f'blast_nt_database = /data/blastdb/nucleotide/nt/{blast_nt.value}/nt',
            f'blast_nr_database = /data/blastdb/protein/nr/{blast_nr.value}/nr'
        ]

    def __build_pfam_options(self):
        """
        Builds the options that specify the location of the pfam database files
        :return: List of options
        """
        pfam = self._parameters['pfam']
        return [
            f'pfamId2TaxId_file = /data/carma/pfam{pfam.value}/pfamid2taxid.txt.gz',
            f'pfamA_txt_file = /data/pfam/{pfam.value}/pfamA_parsed.txt',
            f'gene_ontology_txt_file = /data/pfam/{pfam.value}/parsed_gene_ontology.txt',
            f'pfam_A_hmm_file = /data/pfam/{pfam.value}/Pfam-A.hmm',
            f'pfam_fasta_dir = /data/pfam/{pfam.value}/fasta_dir/'
        ]

    def __build_fixed_options(self):
        """
        Builds the fixed options that specify the location of several binaries and directories
        :return: List of options
        """
        return [f'carma_binary = /usr/local/bin/{self._carma_dir}/carma',
                f'carma_sge_script = /usr/local/bin/{self._carma_dir}/carma3.sh',
                f'hmmfetch_bin = /usr/local/bin/{self._hmmer_dir}/bin/hmmfetch',
                f'hmmalign_bin = /usr/local/bin/{self._hmmer_dir}/bin/hmmalign',
                f'hmmscan_bin = /usr/local/bin/{self._hmmer_dir}/bin/hmmscan',
                'cluster_tmp_dir = /data/temp/',
                'gzip_bin = /bin/gzip',
                'zcat_bin = /bin/zcat']

    @staticmethod
    def __build_unused_options():
        """
        Builds the options that are not implemented in this class as they are experimental but are still required
        in the config file
        :return: List of options
        """
        return [
            'rdp_bact_aligned = /not/implemented',
            'rdp_arch_aligned = /not/implemented',
            'rdp_unaligned = /not/implemented',
            'blast_rdp_database = /not/implemented',
            'blast_rdp_evalue = 0.0'
        ]

    def __get_basename(self):
        """
        Returns the prefix that will be used in the output.
        :return: String with the prefix used in the output
        """
        infile = os.path.basename(self._tool_inputs[self._input_key][0].path)
        return os.path.join(self._folder, os.path.splitext(infile)[0])

    def __get_output_filename(self):
        """
        Returns the output filename based on the input key
        :return: Output filename
        """
        basename = self.__get_basename()
        file_names = {
            'FASTA': basename + '.egt',
            'EGT': basename + '.tax',
            'BLASTX': basename + '.tax'
        }
        return file_names[self._input_key]

    def __get_output_key(self):
        """
        Returns the key for the output file based on the input key
        :return: Output key
        """
        keys = {
            'FASTA': 'EGT',
            'EGT': 'TAX',
            'BLASTX': 'TAX'
        }
        return keys[self._input_key]

    def __set_output(self):
        """
        Sets the name of the output files. If the input key is FASTA, an HMM file will be created, otherwise a taxonomy
        will be created.
        :return: None
        """
        self._tool_outputs[self.__get_output_key()] = [ToolIOFile(Path(self.__get_output_filename()))]

    def __build_input_string(self):
        """
        Creates the string with the input and output files
        :return: String with the input parameters
        """
        items = [f'--input {self._tool_inputs[self._input_key][0].path}',
                 f'--output {self.__get_output_filename()}']
        return ' '.join(items)

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        input_string = self.__build_input_string()
        options_string = ' '.join(self._build_options(self._config_param + self._config_loc_param, ' '))
        options_string += f' --config {self.__create_config_file()}'
        self._command.command = f'{self._tool_command} {input_string} {options_string}'
