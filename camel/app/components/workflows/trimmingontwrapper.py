from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union, Any

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake import snakemakeutils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import trimming_ont
from camel.resources.snakefile.trimming_ont import INPUT_ONT_FASTQ


class TrimmingONTWrapper:
    """
    This class is used as a wrapper class around the ONT read trimming Snakemake workflow.
    """

    @dataclass
    class ReadTrimmingOutput:
        """
        Class that holds the output of the read trimming workflow.
        """
        report_section: HtmlReportSection
        tsv_summary: Path
        trimmed_reads: list[ToolIOFile]
        informs_trimming: dict[str, Any]
        log_file: Optional[Path] = None

    def __init__(self, working_dir: Path) -> None:
        """
        Initializes the read trimming wrapper.
        :param working_dir: Working directory
        """
        self._working_dir = Path(working_dir)
        self._output = None

    def run_workflow(self, se_reads: Path, export_fastq: bool = False, additional_opts: Union[dict, None] = None,
                     threads: int = 8) -> None:
        """
        Runs the read trimming workflow.
        :param se_reads: Input SE FASTQ reads
        :param export_fastq: If True, FASTQ files are included in the report
        :param additional_opts: Additional options
        :param threads: Number of threads to use
        :return: None
        """
        path_io = self._working_dir / INPUT_ONT_FASTQ
        path_io.parent.mkdir(exist_ok=True, parents=True)
        snakemakeutils.dump_object([ToolIOFile(se_reads)], path_io)
        config_data = {
            'working_dir': str(self._working_dir),
            'read_trimming': {
                'export_fastq': str(export_fastq),
                'ont': additional_opts if additional_opts is not None else {}},
            'read_type': 'ont'
        }
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._working_dir)
        output_files = {
            'HTML': trimming_ont.OUTPUT_REPORT,
            'TSV': str(trimming_ont.OUTPUT_SUMMARY).format(ext='tsv'),
            'INFORMS': trimming_ont.OUTPUT_INFORMS,
        }
        SnakePipelineUtils.run_snakemake(
            trimming_ont.SNAKEFILE, config_file, [Path(p) for p in output_files.values()], self._working_dir, threads)
        self.__set_output(output_files)

    def __set_output(self, output_files: dict[str, str]) -> None:
        """
        Sets the output of this tool.
        :param output_files: Output files by key.
        :return: None
        """
        log_path = self._working_dir / 'camel.log'
        self._output = TrimmingONTWrapper.ReadTrimmingOutput(
            report_section=snakemakeutils.load_object(self._working_dir / output_files['HTML'])[0].value,
            tsv_summary=self._working_dir / output_files['TSV'],
            trimmed_reads=snakemakeutils.load_object(
                self._working_dir / trimming_ont.OUTPUT_READS),
            informs_trimming=snakemakeutils.load_object(self._working_dir / output_files['INFORMS']),
            log_file=log_path if log_path.exists() else None
        )

    @property
    def output(self) -> 'TrimmingONTWrapper.ReadTrimmingOutput':
        """
        Returns the report section generated during the read trimming pipeline.
        :return: Trimming section
        """
        return self._output
