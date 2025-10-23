import abc
import argparse
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional, Union

from camel.app.core.errors import SnakemakeExecutionError
from camel.app.core.snakemake import snakepipelineutils
from camel.app.core.utils import fileutils, fastqutils
from camel.app.scriptutils import mainscriptutils, absolute_path_by_pathlib
from camel.app.loggers import fileloggerutils, logger


class BasePipeline(metaclass=abc.ABCMeta):
    """
    Baseclass for pipeline classes.
    """

    def __init__(self, name: str, version: str, snakefile: str, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the pipeline.
        :param name: Pipeline name
        :param version: Pipeline version
        """
        self._name = name
        self._version = version
        self._snakefile = snakefile
        self._args = self._parse_arguments(args)
        self._keep_logs = self._args.log
        self._sample_name = BasePipeline.determine_sample_name(self._args)

    @staticmethod
    @abc.abstractmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments. Should be implemented by the subclasses.
        :param args: Arguments (optional)
        :return: None
        """
        pass

    @staticmethod
    def add_common_arguments(argument_parser: argparse.ArgumentParser) -> None:
        """
        Adds the common arguments for the pipeline
        :param argument_parser: Argument parse
        :return: None
        """
        # Input
        argument_parser.add_argument(
            '--fasta', type=absolute_path_by_pathlib, help="Input FASTA file")
        argument_parser.add_argument(
            '--vcf-unfiltered', type=absolute_path_by_pathlib, help="Input VCF file (only in combination with FASTA file)")
        argument_parser.add_argument('--fasta-name', type=str, help="Input FASTA file name (for Galaxy)")
        argument_parser.add_argument('--sample-name', type=str)
        argument_parser.add_argument(
            '--fastq-pe', nargs=2, type=absolute_path_by_pathlib, help="Input PE FASTQ files")
        argument_parser.add_argument(
            '--fastq-pe-names', nargs=2, help="Input PE FASTQ filenames (for Galaxy)")
        argument_parser.add_argument(
            '--fastq-se', type=absolute_path_by_pathlib, help="Input SE FASTQ file")
        argument_parser.add_argument(
            '--fastq-se-name', help="Input SE FASTQ filename (for Galaxy)")
        argument_parser.add_argument(
            '--input-type', help='Input type',
            choices=['illumina', 'iontorrent', 'ont', 'hybrid', 'fasta', 'fasta_with_vcf'], default='illumina')

        # Output
        argument_parser.add_argument(
            '--working-dir', type=absolute_path_by_pathlib, default=Path(Path.cwd(), 'working'))

        # Options
        argument_parser.add_argument('--threads', default=8, type=int)
        argument_parser.add_argument(
            '--library', help="Adapter library that was used for the sequencing",
            choices=['NexteraPE', 'TruSeq2', 'TruSeq3'], default='NexteraPE')
        argument_parser.add_argument(
            '--trimming-method', help='Trimming method', choices=['trimmomatic', 'fastp'], default='fastp')
        argument_parser.add_argument('--ont-min-qual', default=10, help='Minimum median quality for ONT input data')
        argument_parser.add_argument('--ont-min-length', default=1000, help='Minimum read length for ONT input data')

        # Logging
        argument_parser.add_argument(
            '--galaxy-job-id', type=str, help='Job id of the run in galaxy (used for logging')
        argument_parser.add_argument(
            '--log', action='store_true', help="If this flag is set, config file and error logs are kept")

    @staticmethod
    def determine_sample_name(args: argparse.Namespace) -> str:
        """
        Determines the sample name from the provided arguments
        :return: Sample name
        """
        if args.sample_name is not None:
            return fileutils.make_valid(args.sample_name)
        # FASTA input
        elif args.input_type in ('fasta', 'fasta_with_vcf'):
            if args.fasta_name is not None:
                return fileutils.make_valid(Path(args.fasta_name).stem)
            return fileutils.make_valid(Path(args.fasta).stem)
        # PE reads (illumina / hybrid)
        elif args.input_type in ('illumina', 'hybrid'):
            if args.fastq_pe_names is not None:
                return fastqutils.get_sample_name(args.fastq_pe_names[0], fastqutils.PATTERN_FQ_PE)
            return fastqutils.get_sample_name(args.fastq_pe[0], fastqutils.PATTERN_FQ_PE)
        # SE reads (ont / iontorrent)
        elif args.input_type == 'ont':  # ('ont', 'iontorrent'):
            if args.fastq_se_name is not None:
                return fastqutils.get_sample_name(args.fastq_se_name, fastqutils.PATTERN_FQ_ONT)
            return fastqutils.get_sample_name(args.fastq_se, fastqutils.PATTERN_FQ_ONT)
        raise ValueError('Cannot determine sample name')

    @property
    def name(self) -> str:
        """
        Returns the pipeline name.
        :return: Name
        """
        return self._name

    @property
    def title(self) -> str:
        """
        Returns the title of this pipeline, defaults to pipeline name.
        :return: Title
        """
        return self.name

    @property
    def version(self) -> str:
        """
        Returns the pipeline version.
        :return: Version
        """
        return self._version

    @property
    def sample_name(self) -> str:
        """
        Returns the sample name.
        :return: Sample name
        """
        return self._sample_name

    @property
    def galaxy_job_id(self) -> Union[int, None]:
        """
        Returns the galaxy job id (if there is one).
        :return: Galaxy job id
        """
        return self._args.galaxy_job_id if 'galaxy_job_id' in self._args else None

    def _get_input_links(self) -> list[list[tuple[str, Path, str]]]:
        """
        Returns the links to the input FASTQ files.
        :return: Links (key, path, name)
        """
        links = []

        # FASTA input
        if self._args.input_type in ('fasta', 'fasta_with_vcf'):
            links.append(['fasta', self._args.fasta, f'{self.sample_name}.fasta'])

        # FASTA + VCF input (Mycobacterium HERA)
        if self._args.input_type == 'fasta_with_vcf':
            links.append(['vcf_unfiltered', self._args.vcf_unfiltered, f'{self.sample_name}.vcf'])

        # PE reads
        if self._args.input_type in ('illumina', 'hybrid'):
            for read_nb, path in enumerate(self._args.fastq_pe, start=1):
                gzipped = fileutils.is_gzipped(path)
                links.append(['fastq_pe', path, f"{self.sample_name}_{read_nb}.fastq{'.gz' if gzipped else ''}"])

        # SE reads
        if self._args.input_type in ('hybrid', 'ont'):
            gzipped = fileutils.is_gzipped(self._args.fastq_se)
            links.append(['fastq_se', self._args.fastq_se, f"{self.sample_name}.fastq{'.gz' if gzipped else ''}"])

        # Check if links were created
        if len(links) == 0:
            raise ValueError(f'Invalid input files for input type: {self._args.input_type}')
        return links

    def _symlink_input(self) -> dict[str, list[dict[str, Any]]]:
        """
        Symlinks the input files.
        :return: List of FASTQ input dictionaries
        """
        # Determine link names
        links = self._get_input_links()

        # Create directory
        dir_links = self._args.working_dir / 'input'
        if not dir_links.exists():
            dir_links.mkdir(parents=True)

        # Link files
        dict_input = {}
        for key, path_orig, link_name in links:
            # Add key if missing
            if key not in dict_input:
                dict_input[key] = []

            # Create the symlink
            path_new = dir_links / link_name
            logger.debug(f"Symlinking input file: {path_orig} -> {link_name}")
            if path_new.is_symlink():
                path_new.unlink()
            path_new.symlink_to(path_orig)

            # Store the link
            dict_input[key].append({'name': link_name, 'path': str(path_new)})

        # Return output dictionary
        return dict_input

    def _run_snakemake_main(self, config_file: str) -> None:
        """
        Runs the main snakefile for the pipeline.
        :param config_file: Configuration file
        :return: None
        """
        # Store config file
        if self._keep_logs is True:
            fileloggerutils.store_config_file(Path(config_file), self._name, self.galaxy_job_id)

        # Clear existing Galaxy output files when html output is selected
        if 'output_html' in self._args:
            mainscriptutils.prepare_galaxy_output(Path(self._args.output_dir), Path(self._args.output_html))

        # Path to the logfile
        log_file = self._args.working_dir / 'camel.log'
        try:
            snakepipelineutils.run_snakemake(
                self._snakefile, config_file, [], working_dir=self._args.working_dir, threads=self._args.threads)
            logger.info("Pipeline finished successfully")
        except SnakemakeExecutionError as err:
            if self._keep_logs and log_file.exists():
                log_file = fileloggerutils.store_log_file(log_file, self._name, self.galaxy_job_id, True)
                raise RuntimeError(f"Error executing Snakemake. Check log for more information: {log_file}")
            else:
                raise err

        # Copy log file to output directory if that directory is given
        if log_file.exists() and 'output_dir' in self._args:
            shutil.copyfile(str(log_file), str(Path(self._args.output_dir) / 'camel.log'))

        # Store log file
        if self._keep_logs and log_file.exists():
            fileloggerutils.store_log_file(log_file, self._name, self.galaxy_job_id)

    def get_template_data(self, input_dict: dict[str, list[dict[str, str]]]) -> dict[str, Any]:
        """
        Returns the template data that is common to all pipelines.
        :param input_dict: Dictionary with pipeline input files
        :return: Template data
        """
        template_data: dict[str, Any] = {
            'pipeline': {
                'name': self._name,
                'version': f"{self._version}",
                'title': self.title
            },
            'input': input_dict,
            'input_type': self._args.input_type,
            'sample_name': self.sample_name,
            'working_dir': str(self._args.working_dir),
            'read_trimming': {}
        }
        if self._args.trimming_method is not None:
            template_data['read_trimming']['method'] = self._args.trimming_method
        if self._args.input_type == 'ont':
            template_data['read_trimming']['ont'] = {
                'min_length': self._args.ont_min_length,
                'min_qual': self._args.ont_min_qual,
            }
        return template_data
