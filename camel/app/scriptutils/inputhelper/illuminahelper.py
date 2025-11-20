import dataclasses
from pathlib import Path
from typing import Any

from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.reports.htmlreport import HtmlReport
from camel.app.core.utils import fileutils
from camel.app.loggers import logger
from camel.app.scriptutils import model
from camel.app.scriptutils.basepipe.fastqinput import FastqInput
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.inputhelper.inputhelperbase import (
    InputHelperBase,
    TrimmingOpts,
    AssemblyOpts,
)
from camel.app.wrappers.trimmingilluminawrapper import TrimmingIlluminaWrapper


@dataclasses.dataclass(frozen=True)
class IlluminaTrimmingOpts(TrimmingOpts):
    """
    Options for trimming Illumina reads.
    """
    pass

@dataclasses.dataclass(frozen=True)
class IlluminaAssemblyOpts(AssemblyOpts):
    """
    Options for assembly of Illumina reads.
    """
    assembly_kmers: str | None = dataclasses.field(default=None, metadata= {'help': 'Assembly kmers'})

    def to_dict(self) -> dict[str, Any]:
        """
        Returns the options as a dictionary.
        """
        param_data = {}
        if self.assembly_kmers is not None:
            param_data['kmers'] = self.assembly_kmers
        return param_data


class IlluminaHelper(InputHelperBase[IlluminaTrimmingOpts, IlluminaAssemblyOpts]):
    """
    Helper for Illumina input.
    """

    def __symlink_input(self, fq_in: ScriptInput) -> FastqInput:
        """
        Symlinks the input files.
        :param fq_in: FASTQ input
        :return: Symlinked FASTQ input
        """
        dir_symlinks = self._dir / 'symlinks'
        dir_symlinks.mkdir(exist_ok=True, parents=True)
        paths_symlinks = []
        for path_fq, ori in zip(fq_in.fastq_pe, ('1', '2')):
            path_symlink = dir_symlinks / f"{self._name}_{ori}.fastq{'.gz' if fileutils.is_gzipped(path_fq) else ''}"
            path_symlink.symlink_to(path_fq)
            paths_symlinks.append(path_symlink)
        return FastqInput('illumina', pe=[ToolIOFile(p) for p in paths_symlinks], is_pe=True)

    def prepare_fasta_input(self, script_in: ScriptInput, report: HtmlReport) -> Path:
        """
        Prepares the FASTA input.
        :param script_in: Script input
        :param report: HTML report
        :return: Path to FASTA file
        """
        logger.info("Preparing FASTA input (Illumina data)")

        # ONT FASTQ input
        fq_in_assembly = self.__symlink_input(script_in)
        if self._trimming_opts.trim_reads is True:
            fq_in_assembly = self._trim_reads(fq_in=fq_in_assembly, report=report)

        # De novo assembly
        return self.assemble_fastq_reads(
            fq_in_assembly,
            report=report,
            opts=self._assembly_opts,
            threads=self._threads)

    def _trim_reads(self, fq_in: FastqInput, report: HtmlReport) -> FastqInput:
        """
        Trims the input reads.
        :param fq_in: Script input
        :param report: HTML output report
        :return: FastqInput object
        """
        logger.info("Trimming reads (ONT data)")

        # Run workflow
        trimming = TrimmingIlluminaWrapper(self._dir / 'trimming')
        trimming.run(
            pe_reads=[io.path for io in fq_in.pe],
            method='fastp',
            export_fastq=self._trimming_opts.include_fastq,
            threads=self._threads)

        # Save output
        report_section_trimming = trimming.output.report_section
        report.add_html_object(report_section_trimming)
        report_section_trimming.copy_files(report.output_dir)
        self._logs['trimming'] = trimming.output.log_file
        self._informs.append(trimming.output.informs_trimming)
        return FastqInput('illumina', pe=trimming.output.trimmed_reads_pe, is_pe=True)

    def prepare_fastq_input(self, script_in: ScriptInput, report: HtmlReport) -> FastqInput:
        """
        Prepares the FASTQ input.
        :param script_in: Script input
        :param report: HTML report
        :return: None
        """
        logger.info("Trimming reads (Illumina data)")
        if script_in.type_ is not model.InputType.ILLUMINA:
            raise ValueError(f"Invalid input type: {script_in.type_}")

        # Trim the reads (if needed)
        fq_in = self.__symlink_input(script_in)
        if self._trimming_opts.trim_reads is True:
            return self._trim_reads(fq_in, report=report)
        return fq_in

    @staticmethod
    def opts_from_cli(opts: dict[str, Any]) -> tuple[IlluminaTrimmingOpts, IlluminaAssemblyOpts, int]:
        """
        Parses the options from the command line.
        """
        return (
            IlluminaTrimmingOpts(
                trim_reads=opts['trim_reads'],
                include_fastq=opts['include_fastq'],
            ),
            IlluminaAssemblyOpts(assembly_min_contig_len=opts['assembly_min_contig_len']),
            opts['threads']
        )
