import argparse
import bz2
import datetime
import os
import subprocess
import sys

import yaml

from camel.app.core.command import Command
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.pipeline.pipeline import Pipeline


class GATKSomaticMain:
    """
    Main class to run the GATK somatic variant caller pipeline.
    Generates a config yml file based on CL arguments;
    Dumps a compressed bz2 file of the config file for easy rerun of the pipepline if needed
    (Re)runs the pipeline.
    """
    DB_LOGGING = True
    LOGGING_LEVEL = 'pipeline'
    SNAKEFILE = os.path.join(os.path.dirname(__file__), 'gatk_somatic_steps.snakefile')
    CORES = 5
    TOOL_PARAM_DIR = os.path.join(os.path.dirname(__file__), 'tool_data')

    def __init__(self):
        """
        Initialise.
        :return: None
        """
        self._args = None
        self._ap = None
        self._config_data = dict()
        self._pipeline = None

        self._camel = Camel(tool_parameter_loc=self.TOOL_PARAM_DIR)

        # galaxy dump directory in case of failure
        self._galaxy_dump_dir = os.path.join(self._camel.config["galaxy"]["dump_dir"], "GATK_somatic_calling")

        # dump directory for compressed config files
        self._config_dump_dir = os.path.join(self._camel.config["config_dump_dir"])
        # self._config_dump_dir = "/scratch/todel/temp/GATK_MuTect_configs"

    def run(self) -> None:
        """
        Sets-up and runs the pipeline and logs stderr and stdout if running the command fails.
        Logs the initial input.
        :return: None
        """
        self._pipeline = Pipeline(name='GATK somatic calling', camel=self._camel, logging_level=self.LOGGING_LEVEL)
        self._args = self.__parse_command_line()
        self.__check_inputs()

        if self._args.input_config_file is not None:
            try:
                self.__generate_config_from_previous()
            except FileNotFoundError as error:
                sys.exit("run_gatk_somatic_pipeline.py: error: {}".format(error))
        else:
            self.__generate_config_file()

        self.__archive_config_file()

        self.__log_initial_input()

        self.__execute_snakemake_pipeline()

    def __log_initial_input(self):
        """
        Logs initial inputs.
        :return: 
        """
        input_dict = dict()
        if self.DB_LOGGING:
            if self._args.paired_end:
                input_dict['FASTQ_PE'] = [ToolIOFile(f) for f in self._config_data['fastq']]
            elif self._args.single_end:
                input_dict['FASTQ_SE'] = [ToolIOFile(f) for f in self._config_data['fastq']]

            input_dict['YAML_config'] = [ToolIOFile(self.runtime_config_path)]

            self._pipeline.set_initial_input(input_dict)

    def __execute_snakemake_pipeline(self):
        """
        Execute the snakemake workflow 
        if pipeline is run from galaxy, log stdout and stderr when command fails.
        :return: 
        """
        snakemake_params = " ".join([self._args.unlock, self._args.dag, self._args.dryrun])
        to_execute = 'snakemake --configfile {config} --snakefile {snakefile} --cores {cores} {snakemake_params}'.format(config=self.runtime_config_path, snakefile=self.SNAKEFILE, cores=self._args.threads, snakemake_params=snakemake_params)
        command = Command(to_execute)
        command.run_command(self._args.work_dir, subprocess.STDOUT)
        if command.returncode != 0 and self._args.from_galaxy:
            with open(os.path.join(self._galaxy_dump_dir, "{}.out".format(self._args.job_id)), "w") as file_out:
                file_out.write(command.stdout)
            with open(os.path.join(self._galaxy_dump_dir, "{}.err".format(self._args.job_id)), "w") as file_out:
                file_out.write(command.stderr)
            raise RuntimeError(
                "Error executing Snakemake. Check log ('{}') for more information.".format(os.path.join(self._galaxy_dump_dir, "{}.err".format(
                    self._args.job_id))))

    def __parse_command_line(self):
        """
        Parses the command line arguments.
        :return: argparse.ArgumentParser object
        """
        self._ap = argparse.ArgumentParser(description='Run the GATK somatic variant caller pipeline.')

        # various
        self._ap.add_argument('-W', '--wdir', dest='work_dir', metavar='work_dir', help='Working directory')

        # input
        gp = self._ap.add_mutually_exclusive_group()
        gp.add_argument('-PE', '--paired_end', metavar='fq_file', dest='paired_end', help='Paired-end fastq files.',
                        nargs='+')
        gp.add_argument('-SE', '--single_end', metavar='fq_file', dest='single_end', help='Single-end fastq files.',
                        nargs='+')

        # Variant caller to use
        self._ap.add_argument('-V', '--variant_caller', metavar='variant_caller', dest='variant_caller', help='Variant caller to use.', choices=["mutect1", "mutect2"])


        # output
        self._ap.add_argument('--mutect1_vcf_output', dest='mutect1_vcf_output', metavar='mutect1_vcf_output', help='Output vcf file from MuTect1.')
        self._ap.add_argument('--mutect1_tab_output', dest='mutect1_tab_output', metavar='mutect1_tab_output', help='Output variant-call table file for MuTect1.')
        self._ap.add_argument('--mutect2_vcf_output', dest='mutect2_vcf_output', metavar='mutect2_vcf_output', help='Output vcf file from MuTect2.')
        self._ap.add_argument('--mutect2_bam_output', dest='mutect2_bam_output', metavar='mutect2_bam_output',
                              help='Output bam file from MuTect2: File to which assembled haplotypes should be written.')

        self._ap.add_argument('--covar_output', dest='covar_output', metavar='covar_output', help='Output covariates analysis pdf')
        self._ap.add_argument('--bam_output', dest='bam_output', metavar='bam_output', help='Aligned BAM file')

        # references
        self._ap.add_argument('-R', '--fasta_ref', dest='fasta_ref', metavar='fasta_ref',
                              help='Human genome reference fasta file name (as in db_loc).', choices=["broad_b37_human_Genome_1K_v37"])
        self._ap.add_argument('-S', '--vcf_snps', metavar='vcf_snps', dest='vcf_known_snps',
                              help='Known variant sites (snps) vcf file name (as in db_loc).', choices=["broad_b37_snps_high_confidence"])
        self._ap.add_argument('-I', '--vcf_indels', metavar='vcf_indels', dest='vcf_known_indels',
                              help='Known variant sites (indels) vcf file name (as in db_loc).', choices=["broad_b37_indels_gold_standard"])

        # MarkDuplicates flag
        self._ap.add_argument('--mark_duplicates', dest='markduplicates', help='Mark duplicate reads.', action='store_true')

        # Downsampling
        self._ap.add_argument('--downsampling_type', dest='downsampling_type',
                              help='Type of downsampling to perform on reads (by MuTect). NONE,ALL_READS,BY_SAMPLE. Default: BY_SAMPLE. Perform or not downsampling on reads. '
                             'By default, MuTect downsamples to 1000 reads. Usage example: --downsample None (disables downsampling).')
        self._ap.add_argument('--downsampling_target', dest='downsampling_target', help='Target value for downsampling to perform on reads (by MuTect/MuTect2). Default: 1000. Usage example: --downsample 100000 (sets target value to 100000 reads, effectively disables it).')
        # MuTect1:
        # gap_events_threshold
        self._ap.add_argument('--gap_events_threshold', dest='gap_events_threshold', help='For MuTect1; number of reads allowed to contain insdels around a fixed window (MuTect1 default 11 bp) before being marked as gap_event and filtered-out.')

        # strand_artifact_lod
        self._ap.add_argument('--strand_artifact_lod', dest='strand_artifact_lod', help='For MuTect1; log-odds ratio for strand bias. Default MuTect: 2.0; disable: -99999')

        # MuTect2:
        # run from galaxy flag
        self._ap.add_argument('--from_galaxy', dest='from_galaxy', help='Indicates that the command is run from galaxy. Useful for logging stderr.', action='store_true')

        # number of threads
        self._ap.add_argument('--threads', dest='threads', help='Maximum number of threads for snakemake to allow.', default=self.CORES)

        # downsampling type to perform
        self._ap.add_argument('--mutect2_downsampling_type', dest='MuTect2_downsampling_type',
                        help='Type of downsampling to perform on reads (by MuTect2). NONE,ALL_READS,BY_SAMPLE. Default: NONE.'
                             'Usage example: --downsample None (disables downsampling).', choices=['NONE', 'SAMPLE', 'ALL_READS'])

        # MuTect2 debugging
        # force active regions
        self._ap.add_argument('--mutect2_force_active', dest='MuTect2_force_active',
                        help='Force active regions (see --forceActive param in MuTect2 online doc).', action='store_true')
        self._ap.add_argument('--mutect2_disable_optimizations', dest='MuTect2_disable_optimizations',
                        help='disable optimizations in active regions.', action='store_true')

        # active region output file.
        self._ap.add_argument('--mutect2_output_active_regions', dest='MuTect2_output_active_region_igv',
                        help='Output active region igv file.', action='store_true')
        self._ap.add_argument('--mutect2_active_region_file', dest='MuTect2_active_region_igv_file',
                        help='Output active region igv file name.', default='None')
        # downsampling type to perform
        self._ap.add_argument('--mutect2_output_mode', dest='MuTect2_output_mode',
                        help='Output_mode for vcf (MuTect2). EMIT_VARIANTS_ONLY,EMIT_ALL_CONFIDENT_SITES,EMIT_ALL_SITES.'
                             'Usage example: --downsample None (disables downsampling).', choices=['EMIT_VARIANTS_ONLY', 'EMIT_ALL_CONFIDENT_SITES', 'EMIT_ALL_SITES'])

        # remove first bases in reads
        self._ap.add_argument('--rm_bases_from_reads', dest='bases_to_rm_from_reads',
                              help="Number of bases to systematically remove from reads before procesing. Useful for "
                                   "amplicon-based sequencing if many False Positive technical artifacts are observed at that position. By default, inactive.",
                              default=0)

        # annotation
        self._ap.add_argument('--vep', dest='run_vep', help='Run Variant Effect Predictor; only executed with mutect1.', action='store_true')
        self._ap.add_argument('--vep_output', dest='vep_output', help="Output file name for Variant Effect Predictor.", default="vep_output.vcf")
        self._ap.add_argument('--oncotator', dest='run_oncotator', help='Run oncotator; only executed with mutect1.', action='store_true')
        self._ap.add_argument('--oncotator_output', dest='oncotator_output', help="Output file name for oncotator.", default="oncotator_output.tsv")
        self._ap.add_argument('--oncotator_tx_mode', dest='oncotator_tx_mode', help="tx mode for oncotator (canonical or effect).",
                              default="EFFECT", choices=["CANONICAL", "EFFECT"])

        # snakemake arguments
        # snakemake unlock
        self._ap.add_argument('--unlock', dest='unlock', action="store_const", const="--unlock", help='Unlock snakemake working directory.', default="")
        # snakemake dag
        self._ap.add_argument('--dag', dest='dag', action="store_const", const="--dag", help='Generate snakemake DAG.', default="")
        # dryrun
        self._ap.add_argument('--dryrun', dest='dryrun', action="store_const", const="--dryrun", help='Snakemake dryrun; generates the list of jobs.', default="")

        # job id
        self._ap.add_argument('--job_id', dest='job_id', metavar='job_id', help='Job ID for debugging and logger.',
                              default=datetime.datetime.now().strftime("%Y%m%d_%H%M%S-%f"))
        # mutect threads (mainly for testing)
        self._ap.add_argument('--test_mutect_nct', dest='test_mutect_nct', metavar='test_mutect_nct', help='Threads to use for mutect2. Multithreading leads to longer runtime on targeted data.')

        # re-use pre-generated config file
        self._ap.add_argument('--config_file', dest='input_config_file', metavar='input_config_file', help='Pre-generated config file to use for pipeline run. All other arguments will be not used.', default=None)

        return self._ap.parse_args()

    def __check_inputs(self):
        """
        Checks the inputs provided to the pipeline via the command-line.
        Checks that the inputs are coherent and that required arguments are provided. Extension to argparse.
        :return: 
        """
        required_args_msg = "Required arguments: fasta_ref, vcf_snps, vcf_indels, variant_caller, [paired_end or single_end]."
        try:
            if self._args.input_config_file is None:
                if self._args.fasta_ref is None or self._args.vcf_known_snps is None or self._args.vcf_known_indels is None or self._args.variant_caller is None or (self._args.paired_end is None and self._args.single_end is None):
                    raise ValueError("input_config_file not provided and at least one of the required arguments is missing.")
            else:
                print("Warning: using pre-generated config file.\nThis bypasses the command-line arguments checks and can cause unexpected behaviour. \nCL parameters will be ignored and config file will be used instead.")
        except ValueError as error:
            self._ap.print_usage()
            print("run_gatk_somatic_pipeline.py: error: {}".format(error))
            print(required_args_msg)
            sys.exit(0)

    def __generate_config_file(self):
        """
        Generates a yaml config file based on CLA.
        :return: None
        """

        # Name of config file generated at runtime for snakemake pipeline
        self.runtime_config_path = os.path.join(os.getcwd(), 'runtime_config_{}.yaml'.format(self._args.job_id))

        # logging level
        self._config_data['logging_level'] = self.LOGGING_LEVEL

        # tool parameters directory
        self._config_data['TOOL_PARAM_DIR'] = self.TOOL_PARAM_DIR

        # Add the job id to the config
        self._config_data['pipeline_job_id'] = self._pipeline.job_id
        self._config_data['pipeline_name'] = self._pipeline.name

        # Set working directory: if not given as CLA, set to current working directory (Galaxy)
        if self._args.work_dir:
            self._config_data['working_dir'] = self._args.work_dir
        else:
            self._config_data['working_dir'] = os.getcwd()

        # References
        self._config_data['fasta_ref'] = self._args.fasta_ref
        self._config_data['vcf_known_snps'] = self._args.vcf_known_snps
        self._config_data['vcf_known_indels'] = self._args.vcf_known_indels

        # Input fastq files
        if self._args.paired_end:
            self._config_data['fastq'] = self._args.paired_end
            self._config_data['PE'] = True
            self._config_data['SE'] = False
        if self._args.single_end:
            self._config_data['fastq'] = self._args.single_end
            self._config_data['SE'] = True
            self._config_data['PE'] = False

        # variant caller choice
        self._config_data['variant_caller'] = self._args.variant_caller

        # Output filenames
        if self._args.mutect1_vcf_output:
            self._config_data['mutect1_vcf_output'] = self._args.mutect1_vcf_output
        if self._args.mutect1_tab_output:
            self._config_data['mutect1_tab_output'] = self._args.mutect1_tab_output
        if self._args.mutect2_vcf_output:
            self._config_data['mutect2_vcf_output'] = self._args.mutect2_vcf_output
        if self._args.mutect2_bam_output:
            self._config_data['mutect2_bam_output'] = self._args.mutect2_bam_output
        if self._args.covar_output:
            self._config_data['covar_output'] = self._args.covar_output
        if self._args.bam_output:
            self._config_data['bam_output'] = self._args.bam_output

        # Flag for MarkDuplicates
        self._config_data['run_markDuplicates'] = self._args.markduplicates

        # MuTect1 parameters
        if self._config_data['variant_caller'] == "mutect1":
            # Downsampling
            if self._args.downsampling_type:
                self._config_data['MuTect1_downsampling_type'] = self._args.downsampling_type
            if self._args.downsampling_target:
                self._config_data['MuTect1_downsampling_target'] = self._args.downsampling_target
            # gap_event_threshold
            if self._args.gap_events_threshold:
                self._config_data['gap_events_threshold'] = self._args.gap_events_threshold
            # strand_artifact_lod
            if self._args.strand_artifact_lod:
                self._config_data['strand_artifact_lod'] = self._args.strand_artifact_lod

        # MuTect2 parameters
        if self._config_data['variant_caller'] == "mutect2":
            # Downsampling
            if self._args.downsampling_target:
                self._config_data['MuTect2_downsampling_target'] = self._args.downsampling_target
            if self._args.MuTect2_downsampling_type:
                self._config_data['MuTect2_downsampling_type'] = self._args.MuTect2_downsampling_type
            if self._args.MuTect2_force_active:
                self._config_data['MuTect2_force_active'] = True
            if self._args.MuTect2_disable_optimizations:
                self._config_data['MuTect2_disable_optimizations'] = True
            if self._args.MuTect2_output_active_region_igv:
                self._config_data['MuTect2_output_active_region_igv'] = self._args.MuTect2_output_active_region_igv
            if self._args.MuTect2_active_region_igv_file:
                self._config_data['MuTect2_active_region_igv_file'] = self._args.MuTect2_active_region_igv_file

            if self._args.MuTect2_output_mode:
                self._config_data['MuTect2_output_mode'] = self._args.MuTect2_output_mode

        # Indel realigner flag (indel realignment only required for MuTect1)
        if self._args.variant_caller == "mutect1":
            self._config_data['run_indel_realignment'] = True
        elif self._args.variant_caller == "mutect2":
            self._config_data['run_indel_realignment'] = False

        # test mutect threads
        if self._args.test_mutect_nct:
            self._config_data['mutect_nct'] = self._args.test_mutect_nct

        # threads
        if self._args.threads:
            self._config_data['threads'] = self._args.threads
        else:
            self._config_data['threads'] = self._args.threads

        # run from Galaxy (for output files naming)
        self._config_data['from_galaxy'] = self._args.from_galaxy

        # remove bases from beginning of reads
        if self._args.bases_to_rm_from_reads == "0":
            self._config_data['run_remove_first_base'] = False
        else:
            self._config_data['run_remove_first_base'] = True
            self._config_data['bases_to_rm_from_reads'] = self._args.bases_to_rm_from_reads

        # annotations (Beta; only executed with mutect1 for now)
        if self._args.variant_caller == "mutect1":
            if self._args.run_vep:
                self._config_data['run_vep'] = True
                self._config_data['out_vep'] = self._args.vep_output
            else:
                self._config_data['run_vep'] = False

            if self._args.run_oncotator:
                self._config_data['run_oncotator'] = True
                self._config_data['out_oncotator'] = self._args.oncotator_output
                self._config_data['oncotator_tx_mode'] = self._args.oncotator_tx_mode
            else:
                self._config_data['run_oncotator'] = False
        else:
            if self._args.run_oncotator:
                print("--oncotator parameter set, but ignored because not using mutect1.")
            if self._args.run_vep:
                print("--vep parameter set, but ignored because not using mutect1.")
            self._config_data['run_vep'] = False
            self._config_data['run_oncotator'] = False

        # Create and write to config file
        with open(self.runtime_config_path, 'w') as handle:
            yaml.dump(self._config_data, handle)

    def __generate_config_from_previous(self):
        """
        Imports data from pre-gerenerated config yaml file (replaces previous data) and 
        generates new config file with new pipeline_job_id and runtime_config filename.
        If input file doesn't exist, raise error.
        :return: 
        """
        previous_config_path = os.path.join(os.getcwd(), self._args.input_config_file)
        if os.path.isfile(previous_config_path):
            with open(previous_config_path, "r") as handle_in:
                self._config_data = yaml.load(handle_in)
        else:
            raise FileNotFoundError("Config file '{}' not found.".format(previous_config_path))

        # pipeline job id
        self._config_data['pipeline_job_id'] = self._pipeline.job_id

        # Name of config file generated at runtime for snakemake pipeline
        if self._args.job_id is not None:
            self.runtime_config_path = os.path.join(os.getcwd(), 'runtime_config_{}.yaml'.format(self._args.job_id))
        else:
            self.runtime_config_path = os.path.join(os.getcwd(),
                                                    'runtime_config_{}.yaml'.format(datetime.datetime.now().strftime("%Y%m%d_%H%M%S-%f")))

        # Create and write to config file
        with open(self.runtime_config_path, 'w') as handle:
            yaml.dump(self._config_data, handle)

    def __archive_config_file(self):
        """
        Dump yaml config file as compressed bz2 file in the config dump dir.
        :return: 
        """
        outfile = os.path.join(self._config_dump_dir, 'config_{}_{}.bz2'.format(self._args.job_id, datetime.datetime.now().strftime("%Y%m%d_%H%M%S-%f")))
        with bz2.BZ2File(outfile, 'wb', compresslevel=9) as handle:
            yaml.dump(self._config_data, handle, encoding=('utf-8'))


if __name__ == '__main__':
    main = GATKSomaticMain()
    main.run()
