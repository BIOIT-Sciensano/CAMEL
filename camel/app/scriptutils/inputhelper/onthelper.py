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
    AssemblyOpts
)
from camel.app.tools.seqtk.seqtkconvert import SeqtkConvert
from camel.app.wrappers.trimmingontwrapper import TrimmingONTWrapper


@dataclasses.dataclass(frozen=True)
class ONTTrimmingOpts(TrimmingOpts):
    """
    Options for trimming ONT reads.
    """
    ont_min_len: int | None = dataclasses.field(default=500, metadata={
        'help': 'Minimum length of ONT reads',
        'show_default': True,
    })
    ont_min_qual: int | None = dataclasses.field(default=10, metadata={
        'help': 'Minimum quality of ONT reads',
        'show_default': True
    })

@dataclasses.dataclass(frozen=True)
class ONTAssemblyOpts(AssemblyOpts):
    """
    Options for assembly of ONT reads.
    """
    assembly_flye_meta: bool = False

    def to_dict(self) -> dict[str, Any]:
        """
        Returns the options as a dictionary.
        """
        return {'meta': self.assembly_flye_meta}


class ONTHelper(InputHelperBase[ONTTrimmingOpts, ONTAssemblyOpts]):
    """
    Helper for ONT input.
    """

    def __symlink_input(self, script_in: ScriptInput) -> FastqInput:
        """
        Symlinks the input files.
        :param script_in: Script input
        :return: Symlinked FASTQ input
        """
        dir_symlinks = self._dir / 'symlinks'
        dir_symlinks.mkdir(exist_ok=True, parents=True)
        path_symlink = dir_symlinks / f"{self._name}.fastq{'.gz' if fileutils.is_gzipped(script_in.fastq_se) else ''}"
        path_symlink.symlink_to(script_in.fastq_se)
        return FastqInput('ont', se=[ToolIOFile(path_symlink)], is_pe=False)

    def prepare_fasta_input(self, script_in: ScriptInput, report: HtmlReport) -> Path:
        """
        Prepares the FASTA input.
        :param script_in: Script input
        :param report: HTML report
        :return: Path to FASTA file
        """
        logger.info("Preparing FASTA input (ONT data)")
        if script_in.type_ is not model.InputType.ONT:
            raise ValueError(f"Invalid input type: {script_in.type_}")

        # ONT FASTQ input
        fq_in_assembly = self.__symlink_input(script_in)
        if self._trimming_opts.trim_reads is True:
            fq_in_assembly = self._trim_reads(fastq_in=fq_in_assembly, report=report)

        # De novo assembly
        return self.assemble_fastq_reads(
            fq_in_assembly,
            report=report,
            opts=self._assembly_opts,
            threads=self._threads)

    def _trim_reads(self, fastq_in: FastqInput, report: HtmlReport) -> FastqInput:
        """
        Trims the input reads.
        :param fastq_in: FASTQ input
        :param report: HTML output report
        :return: FastqInput object
        """
        # Run workflow
        trimming = TrimmingONTWrapper(self._dir / 'trimming')
        additional_opts = {}
        if self._trimming_opts.ont_min_len is not None:
            additional_opts['min_length'] = self._trimming_opts.ont_min_len
        if self._trimming_opts.ont_min_qual is not None:
            additional_opts['min_qual'] = self._trimming_opts.ont_min_qual
        trimming.run(
            se_reads=fastq_in.se[0].path,
            export_fastq=self._trimming_opts.include_fastq,
            additional_opts=additional_opts,
            threads=self._threads)

        # Save output
        report_section_trimming = trimming.output.report_section
        report.add_html_object(report_section_trimming)
        report_section_trimming.copy_files(report.output_dir)
        self._logs['trimming'] = trimming.output.log_file
        self._informs.append(trimming.output.informs_trimming)
        return FastqInput('ont', se=trimming.output.trimmed_reads, is_pe=False)

    def prepare_fastq_input(self, script_in: ScriptInput, report: HtmlReport) -> FastqInput:
        """
        Prepares the FASTQ input.
        :param script_in: Script input
        :param report: HTML report
        :return: None
        """
        logger.info("Trimming reads (ONT data)")
        if script_in.type_ is not model.InputType.ONT:
            raise ValueError(f"Invalid input type: {script_in.type_}")

        # Trim the reads (if needed)
        fastq_in = self.__symlink_input(script_in)
        if self._trimming_opts.trim_reads is True:
            return self._trim_reads(fastq_in, report=report)
        return FastqInput('ont', se=fastq_in.se, is_pe=False)

    @staticmethod
    def opts_from_cli(opts: dict[str, Any]) -> tuple[ONTTrimmingOpts, ONTAssemblyOpts, int]:
        """
        Parses the options from the command line.
        :param opts: Command line options
        :return: Assembly options, trimming options, number of threads
        """
        return (
            ONTTrimmingOpts(
                trim_reads=opts['trim_reads'],
                include_fastq=opts['include_fastq'],
                ont_min_len=opts['ont_min_len'],
                ont_min_qual=opts['ont_min_qual'],
            ),
            ONTAssemblyOpts(opts['assembly_flye_meta']),
            opts['threads']
        )

    def prepare_fastq_read_input(self, script_in: ScriptInput, report: HtmlReport) -> Path:
        """
        Prepares the FASTA input by converting FASTQ to FASTA with Seqtk.
        :param script_in: Script input
        :param report: HTML report
        :return: Path to FASTA file
        """
        fastq_in = self.__symlink_input(script_in)
        if self._trimming_opts.trim_reads:
            fastq_in = self._trim_reads(fastq_in, report=report)

        # Convert to FASTA format
        convert = SeqtkConvert()
        convert.add_input_files({'FASTQ': fastq_in.se})
        working_dir = Path(self._dir / 'fq_to_fa')
        working_dir.mkdir(parents=True, exist_ok=True)
        convert.run(working_dir)
        return Path(convert.tool_outputs['FASTA'][0].path)
