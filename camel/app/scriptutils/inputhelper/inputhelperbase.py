import abc
import dataclasses
import shutil
from pathlib import Path
from typing import Any, TypeVar, Generic

from camel.app.core.reports import reportutils
from camel.app.core.reports.htmlreport import HtmlReport
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.loggers import logger
from camel.app.scriptutils.basepipe.fastqinput import FastqInput
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.wrappers.assemblywrapper import AssemblyWrapper


@dataclasses.dataclass(frozen=True)
class TrimmingOpts:
    """
    Dataclass to store trimming options.
    """
    trim_reads: bool = False
    include_fastq: bool = False

@dataclasses.dataclass(frozen=True)
class AssemblyOpts(metaclass=abc.ABCMeta):
    """
    Dataclass to store assembly options.
    """
    assembly_min_contig_len: int = dataclasses.field(default=1000, metadata={'help': 'Minimum contig length for assembly'})

    @abc.abstractmethod
    def to_dict(self):
        """
        Returns the assembly options as a dictionary.
        :return: Options as dictionary
        """
        pass


TTrimmingOpts = TypeVar("TTrimmingOpts", bound=TrimmingOpts)
TAssemblyOpts = TypeVar("TAssemblyOpts", bound=AssemblyOpts)


class InputHelperBase(Generic[TTrimmingOpts, TAssemblyOpts], metaclass=abc.ABCMeta):
    """
    Baseclass for input type helper.
    """

    def __init__(self, dir_: Path, name: str) -> None:
        """
        Initializes the read type helper.
        :param dir_: Working directory
        :param name: Sample name
        :return: None
        """
        self._dir = dir_
        self._name = name
        self._informs: list[dict] = []
        self._logs: dict[str, Path] = {}
        self._trimming_opts: TTrimmingOpts | None = None
        self._assembly_opts: TAssemblyOpts | None = None
        self._threads: int = 1

    @property
    def working_dir(self) -> Path:
        """
        Returns the working directory for the helper.
        :return: Working directory
        """
        return self._dir

    @property
    def logs(self) -> dict[str, Path]:
        """
        Returns the log files (key: name, value: log file path).
        :return: Logs
        """
        return self._logs

    @property
    def informs(self) -> list[dict[str, str]]:
        """
        Returns the informs.
        :return: List of informs
        """
        return self._informs

    def prepare_fasta_input(self, script_in: ScriptInput, report: HtmlReport) -> Path:
        """
        Prepares the FASTA input.
        :param script_in: Script input
        :param report: HTML report
        :return: Path to FASTA file
        """
        pass

    def prepare_fastq_input(self, script_in: ScriptInput, report: HtmlReport) -> FastqInput:
        """
        Prepares the FASTQ input.
        :param script_in: Script input
        :param report: HTML report
        :return: None
        """
        pass

    def set_opts(self, trimming_opts: TTrimmingOpts | None, assembly_opts: TAssemblyOpts | None, threads: int) -> None:
        """
        Sets the helper options.
        :param trimming_opts: Trimming options
        :param assembly_opts: Assembly options
        :param threads: Number of threads
        :return: None
        """
        self._trimming_opts = trimming_opts
        self._assembly_opts = assembly_opts
        self._threads = threads

    def assemble_fastq_reads(
            self, assembly_input: FastqInput, opts: AssemblyOpts, report: HtmlReport | None = None, threads: int = 4) -> Path:
        """
        Assembles FASTQ reads using SPAdes
        :param assembly_input: Assembly input
        :param opts: Assembly options
        :param report: If set, the output is added to the given report
        :param threads: Number of threads
        :return: ToolIOFile FASTA object with the assembled contigs
        """
        logger.info("Starting de novo assembly")
        assembly = AssemblyWrapper(self._dir / 'assembly', assembly_input.read_type)

        # Perform the assembly
        assembly.run(
            self._name,
            fastq_in=assembly_input,
            min_ctg_len=opts.assembly_min_contig_len,
            assembler_opts=opts.to_dict(),
            threads=threads)

        # Save output to the report
        if report is not None:
            report.add_html_object(assembly.output.report_section)
            assembly.output.report_section.copy_files(report.output_dir)
            report.save()

        # Save log file and informs
        if assembly.output.log_file is not None:
            self._logs['assembly'] = assembly.output.log_file
        self._informs.extend(assembly.output.informs)
        return assembly.output.fasta_contigs

    @staticmethod
    @abc.abstractmethod
    def opts_from_cli(opts: dict[str, Any]) -> tuple[TTrimmingOpts | None, TAssemblyOpts | None, int]:
        """
        Parses the options from the command line.
        :param opts: Command line options
        :return: Trimming options, assembly options, number of threads
        """
        raise NotImplementedError()

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
            section_commands = reportutils.create_commands_section(self._informs, self._dir)
            report.add_html_object(section_commands)
        report.save()
