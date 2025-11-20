from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.snakemake import snakemakeutils
from camel.app.core.snakemake import snakepipelineutils
from camel.app.loggers import logger
from camel.snakefiles import trimming_illumina


class TrimmingIlluminaWrapper:
    """
    This class is used as a wrapper class around the Illumina read trimming Snakemake workflow.
    """

    @dataclass
    class ReadTrimmingOutput:
        """
        Class that holds the output of the read trimming workflow.
        """
        report_section: HtmlReportSection
        tsv_summary: Path
        trimmed_reads_pe: list[ToolIOFile]
        trimmed_reads_se_fwd: list[ToolIOFile]
        trimmed_reads_se_rev: list[ToolIOFile]
        informs_trimming: dict[str, Any]
        fastq_reports_pre: list[ToolIOFile]
        fastq_reports_post: list[ToolIOFile]
        log_file: Optional[Path] = None

    def __init__(self, working_dir: Path) -> None:
        """
        Initializes the read trimming wrapper.
        :param working_dir: Working directory
        """
        self._working_dir = Path(working_dir)
        self._output = None

    def run(self, pe_reads: list[Path], adapter: Optional[str] = None, threads: int = 8, export_fastq: bool = False,
            method: Optional[str] = None) -> None:
        """
        Runs the read trimming workflow.
        :param pe_reads: Input PE FASTQ reads
        :param adapter: Adapter to trim
        :param threads: Number of threads to use
        :param export_fastq: If True, FASTQ files are included in the report
        :param method: Trimming method ('trimmomatic' or 'fastp')
        :return: None
        """
        logger.info("Running read trimming workflow (Illumina)")

        # Create the config file
        config_data: dict[str, Any] = {
            'working_dir': str(self._working_dir),
            'input': {'fastq_pe': [{'name': p.name, 'path': str(p)} for p in pe_reads]},
            'read_trimming': {'export_fastq': str(export_fastq)}
        }
        # Add the method (if specified)
        if method is not None:
            config_data['read_trimming']['method'] = method

        if adapter is not None:
            config_data['read_trimming']['adapter'] = adapter
        config_file = snakepipelineutils.generate_config_file(config_data, self._working_dir)

        # Dump the input files in an IO file
        io_pickle_in = Path(self._working_dir / trimming_illumina.INPUT_FASTQ)
        io_pickle_in.parent.mkdir(exist_ok=True, parents=True)
        snakemakeutils.dump_object([ToolIOFile(x) for x in pe_reads], io_pickle_in)

        # Output files
        output_files = {
            'FASTQ': trimming_illumina.OUTPUT_DICT,
            'HTML': trimming_illumina.OUTPUT_REPORT,
            'INFORMS': trimming_illumina.OUTPUT_INFORMS,
            'TSV': str(trimming_illumina.OUTPUT_SUMMARY).format(ext='tsv'),
            'FASTQC_PRE': trimming_illumina.OUTPUT_FASTQC_HTML_PRE,
            'FASTQC_POST': trimming_illumina.OUTPUT_FASTQC_HTML_POST,
        }
        snakepipelineutils.run_snakemake(
            snakefile=trimming_illumina.SNAKEFILE,
            config_path=config_file,
            targets=[Path(p) for p in output_files.values()],
            working_dir=self._working_dir,
            threads=threads)
        self.__set_output(output_files)

    def __set_output(self, output_files: dict[str, str]) -> None:
        """
        Sets the output of this tool.
        :param output_files: Output files by key.
        :return: None
        """
        log_path = self._working_dir / 'camel.log'
        fq_dict = snakemakeutils.load_object(self._working_dir / output_files['FASTQ'])
        self._output = TrimmingIlluminaWrapper.ReadTrimmingOutput(
            report_section=snakemakeutils.load_object(self._working_dir / output_files['HTML'])[0].value,
            tsv_summary=self._working_dir / output_files['TSV'],
            trimmed_reads_pe=fq_dict['PE'],
            trimmed_reads_se_fwd=fq_dict['SE_FWD'] if 'SE_FWD' in fq_dict else [],
            trimmed_reads_se_rev=fq_dict['SE_REV'] if 'SE_REV' in fq_dict else [],
            informs_trimming=snakemakeutils.load_object(self._working_dir / output_files['INFORMS']),
            fastq_reports_pre=snakemakeutils.load_object(self._working_dir / output_files['FASTQC_PRE']),
            fastq_reports_post=snakemakeutils.load_object(self._working_dir / output_files['FASTQC_POST']),
            log_file=log_path if log_path.exists() else None
        )

    @property
    def output(self) -> 'TrimmingIlluminaWrapper.ReadTrimmingOutput':
        """
        Returns the report section generated during the read trimming pipeline.
        :return: Trimming section
        """
        return self._output
