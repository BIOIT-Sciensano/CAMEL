import argparse
from pathlib import Path
from typing import Union

from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.workflows.inputtype.inputtypehelperbase import InputTypeHelperBase
from camel.app.components.workflows.trimmingontwrapper import TrimmingONTWrapper
from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.seqtk.seqtkconvert import SeqtkConvert

class ONTHelper(InputTypeHelperBase):
    """
    Helper class for ONT input.
    """

    def __symlink_ont_reads(self, fastq_file: Union[Path, None], sample_name: str) -> Path:
        """
        Symlinks the input files to a standardized format based on the sample name.
        :param fastq_file: Input ONT FASTQ file
        :param sample_name: Sample name
        :return: Path to renamed file
        """
        if fastq_file is None:
            raise ValueError("ONT data should be SE")
        new_name = f"{sample_name}.fastq{'.gz' if FileSystemHelper.is_gzipped(fastq_file) else ''}"
        return self.symlink_input_files([Path(fastq_file)], [new_name])[0]

    def trim_reads(self, fastq_input: FastqInput, report: HtmlReport, include_fastq: bool, threads: int = 4,
                   **kwargs) -> FastqInput:
        """
        Filters the input ONT reads.
        :param fastq_input: FASTQ input
        :param report: HTML report
        :param include_fastq: Boolean to indicate if FASTQ files should be included in the report
        :param threads: Nb. of threads
        :return: FastqInput object with filtered reads
        """
        logger.info("Trimming reads (ONT data)")
        if fastq_input.is_pe or fastq_input.se is None:
            raise ValueError("ONT input should be SE")

        # Run workflow
        trimming = TrimmingONTWrapper(self._working_dir / 'trimming')
        additional_opts = {}
        if 'min_length' in kwargs:
            additional_opts['min_length'] = kwargs['min_length']
        if 'min_qual' in kwargs:
            additional_opts['min_qual'] = kwargs['min_qual']
        trimming.run_workflow(
            se_reads=Path(fastq_input.se[0].path),
            export_fastq=include_fastq,
            additional_opts=additional_opts,
            threads=threads)

        # Save output
        report_section_trimming = trimming.output.report_section
        report.add_html_object(report_section_trimming)
        report_section_trimming.copy_files(report.output_dir)
        self._log_files['trimming'] = trimming.output.log_file
        self._informs.append(trimming.output.informs_trimming)
        return FastqInput('ont', se=trimming.output.trimmed_reads, is_pe=False)

    def prepare_fasta_input(self, report: HtmlReport, args: argparse.Namespace) -> Path:
        """
        Prepares the FASTA input.
        :param report: HTML report
        :param args: Command-line arguments
        :return: FASTA file
        """
        logger.info("Preparing FASTA input (ONT data)")

        # ONT FASTQ input
        fq_input_se = self.__symlink_ont_reads(args.fastq_se, self._sample_name)
        fastq_input = FastqInput('ont', se=[ToolIOFile(fq_input_se)], is_pe=False)
        if args.trim_reads:
            fastq_input = self.trim_reads(
                fastq_input=fastq_input,
                report=report,
                include_fastq=args.report_include_fastq,
                min_length=args.ont_min_length,
                min_qual=args.ont_min_qual,
                threads=args.threads
            )

        # De novo assembly
        return self.assemble_fastq_reads(fastq_input, args, report)

    def prepare_fastq_input(self, report: HtmlReport, args: argparse.Namespace) -> FastqInput:
        """
        Prepares the FASTQ input.
        :param report: HTML report
        :param args: Command-line arguments
        :return: FASTQ input
        """
        fq_input_se = self.__symlink_ont_reads(Path(args.fastq_se), self._sample_name)
        fastq_input = FastqInput('ont', se=[ToolIOFile(Path(fq_input_se))], is_pe=False)
        if args.trim_reads:
            return self.trim_reads(
                fastq_input=fastq_input,
                report=report,
                include_fastq=args.report_include_fastq,
                min_length=args.ont_min_length,
                min_qual=args.ont_min_qual,
                threads=args.threads
            )
        return FastqInput('ont', se=[ToolIOFile(Path(fq_input_se))], is_pe=False)

    def prepare_fasta_read_input(self, report: HtmlReport, args: argparse.Namespace) -> Path:
        """
        Prepares the FASTA input by converting FASTQ to FASTA with Seqtk.
        :param report: HTML report
        :param args: Command-line arguments
        :return: Path to FASTA file
        """
        fq_input_se = self.__symlink_ont_reads(Path(args.fastq_se), self._sample_name)
        fastq_input = FastqInput(args.input_type, se=[ToolIOFile(Path(fq_input_se))], is_pe=False)
        if args.trim_reads:
            convert_input = self.trim_reads(fastq_input, report, args.report_include_fastq, args.threads)
        else:
            convert_input = fastq_input
        convert = SeqtkConvert()
        convert.add_input_files({'FASTQ': convert_input.se})
        working_dir = Path(self._working_dir / 'conversion')
        working_dir.mkdir(parents=True, exist_ok=True)
        convert.run(working_dir)
        return Path(convert.tool_outputs['FASTA'][0].path)
