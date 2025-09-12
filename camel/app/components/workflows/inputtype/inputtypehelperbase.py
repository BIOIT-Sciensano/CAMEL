import abc
import argparse
import shutil
from pathlib import Path
from typing import Optional, Any

from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.workflows.assemblywrapper import AssemblyWrapper
from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.loggers import logger
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils


class InputTypeHelperBase(metaclass=abc.ABCMeta):
    """
    Base class for input type-specific helper classes.
    """

    def __init__(self, working_dir: Path, sample_name: str) -> None:
        """
        Initializes the helper class.
        :param working_dir: Working directory
        """
        self._working_dir = working_dir
        self._working_dir.mkdir(exist_ok=True, parents=True)
        self._sample_name = sample_name
        self._log_files = {}
        self._informs = []

    @property
    def logs(self) -> dict[str, str]:
        """
        Returns the log files (key: name, value: log file path).
        :return: Logs
        """
        return self._log_files

    @property
    def informs(self) -> list[dict[str, str]]:
        """
        Returns the informs.
        :return: List of informs
        """
        return self._informs

    @property
    def working_dir(self) -> Path:
        """
        Returns the working directory.
        :return: Working directory
        """
        return self._working_dir

    def symlink_input_files(self, files: list[Path], names: Optional[list[str]] = None) -> list[Path]:
        """
        Creates symbolic links for list of files to the working directory.
        :param files: List of input files
        :param names: Target file names
        :return: List of symlink locations
        """
        logger.info(f"Symlinking input files: {', '.join([f.name for f in files])}")
        dir_links = self._working_dir / 'input'
        dir_links.mkdir(exist_ok=True, parents=True)

        # Determine the names
        if names is None:
            names = [None] * len(files)

        linked_files = []
        for file_, name_ in zip(files, names):
            new_symlink = dir_links / FileSystemHelper.make_valid(name_ if name_ is not None else file_.name)
            if new_symlink.exists():
                new_symlink.unlink()
            logger.debug(f"Creating link: {new_symlink.name} -> {file_.name}")
            new_symlink.symlink_to(file_)
            linked_files.append(new_symlink)
        return linked_files

    def __get_assembly_illumina_args(self, args: argparse.Namespace) -> dict[str, Any]:
        """
        Adds the options for the Illumina assembly.
        :param args: Command-line arguments
        :return: None
        """
        options = {}

        # Coverage cut-off
        if args.assembly_cov_cutoff is None:
            options['cov_cutoff'] = 'off'
        elif args.assembly_cov_cutoff == 0:
            options['cov_cutoff'] = 'auto'
        else:
            options['cov_cutoff'] = str(args.assembly_cov_cutoff)

        # Kmers
        if args.assembly_kmers is not None:
            options['kmers'] = args.assembly_kmers
        return options

    def __get_assembly_ont_args(self, args: argparse.Namespace) -> dict[str, Any]:
        """
        Adds the options for the ONT assembly.
        :param args: Command-line arguments
        :return: None
        """
        if args.assembly_flye_meta is True:
            return {'meta': True}
        return {}

    def assemble_fastq_reads(
            self, assembly_input: FastqInput, args: argparse.Namespace, report: Optional[HtmlReport] = None) -> Path:
        """
        Assembles FASTQ reads using SPAdes
        :param assembly_input: Assembly input
        :param report: If set, the output is added to the given report
        :param args: Command line arguments
        :return: ToolIOFile FASTA object with the assembled contigs
        """
        logger.info("Starting de-novo assembly")
        assembly = AssemblyWrapper(self._working_dir / 'assembly', assembly_input.read_type)

        # Additional assembler options
        if assembly_input.read_type == 'illumina':
            assembler_opts = self.__get_assembly_illumina_args(args)
        elif assembly_input.read_type == 'ont':
            assembler_opts = self.__get_assembly_ont_args(args)
        else:
            assembler_opts = {}

        # Perform the assembly
        assembly.run(
            self._sample_name,
            fastq_in=assembly_input,
            min_ctg_len=args.assembly_min_contig_length,
            assembler_opts=assembler_opts,
            threads=args.threads)

        # Save output to the report
        if report is not None:
            report.add_html_object(assembly.output.report_section)
            assembly.output.report_section.copy_files(report.output_dir)
            report.save()

        # Save log file and informs
        if assembly.output.log_file is not None:
            self._log_files['assembly'] = assembly.output.log_file
        self._informs.extend(assembly.output.informs)
        return assembly.output.fasta_contigs

    def export_log_files(self, output_dir: Path) -> None:
        """
        Exports the log files to the output directory.
        :param output_dir: Output directory
        :return: None
        """
        dir_logs = output_dir / 'logs'
        dir_logs.mkdir(parents=True, exist_ok=True)
        for key, path in self.logs.items():
            if path is None:
                logger.warning(f'No log file found for: {key}')
                continue
            shutil.copyfile(path, str(dir_logs / f'log_{key}.txt'))

    def export_output_and_commands_section(self, report: HtmlReport, section: HtmlReportSection) -> None:
        """
        Adds the output and commands sections to the report.
        Copies the log files to the output folder.
        :param report: Report
        :param section: Section to add
        :return: None
        """
        report.add_html_object(section)
        section.copy_files(report.output_dir)
        self.export_log_files(Path(report.output_dir))
        if len(self._informs) > 0:
            section_commands = SnakePipelineUtils.create_commands_section(self._informs, self._working_dir)
            report.add_html_object(section_commands)
        report.save()

    @abc.abstractmethod
    def trim_reads(self, fastq_input: FastqInput, report: HtmlReport, include_fastq: bool, threads: int, **kwargs) -> \
            FastqInput:
        """
        Base function for read-type specific trimming.
        :return: Trimmed FASTQ files
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def prepare_fasta_input(self, report: HtmlReport, args: argparse.Namespace) -> Path:
        """
        Prepares the FASTA input.
        :param report: HTML report
        :param args: Command-line arguments
        :return: FASTA file
        """
        raise NotImplementedError()
