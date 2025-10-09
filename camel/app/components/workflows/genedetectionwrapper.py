from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

from camel.app.components.genedetection.genedetectionhitbase import GeneDetectionHitBase
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake import snakemakeutils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import gene_detection


@dataclass
class GeneDetectionOutput:
    """
    Holder for the gene detection output.
    """
    report_section: HtmlReportSection
    detected_hits: list[GeneDetectionHitBase]
    informs: dict[str, Any]
    log_file: Optional[Path] = None


class GeneDetectionWrapper:
    """
    This class is used as a wrapper class around the gene detection Snakemake workflow.
    """

    def __init__(self, working_dir: Path) -> None:
        """
        Initializes the read trimming helper.
        :param working_dir: Working directory
        """
        self._working_dir = working_dir
        self._output = None

    def __run_workflow(self, config_data: dict[str, Any], threads: int) -> None:
        """
        Runs the gene detection workflow with the given config data and number of threads.
        :param config_data: Config data
        :param threads: Number of threads
        :return: None
        """
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._working_dir)
        output_files = {
            'report': gene_detection.OUTPUT_REPORT.format(db='db'),
            'informs': gene_detection.OUTPUT_INFORMS.format(db='db'),
            'hits': gene_detection.OUTPUT_ALL_HITS.format(db='db')
        }
        SnakePipelineUtils.run_snakemake(
            gene_detection.SNAKEFILE, config_file, [Path(x) for x in output_files.values()], self._working_dir, threads)
        self.__set_output(output_files)

    def run_workflow_blast(self, fasta_path: Path, sample_name: str, db_data: dict[str, Any], threads: int = 8) -> None:
        """
        Runs the gene detection workflow using BLAST.
        :param fasta_path: Input FASTA file
        :param sample_name: Sample name
        :param db_data: Database configuration data (should contain at least 'path' as key referring to the location of
          the database.)
        :param threads: Number of threads to use
        :return: None
        """
        self.__create_input_blast(fasta_path)
        config_data = self.__get_config_data(sample_name, db_data, 'blast', 'illumina')
        self.__run_workflow(config_data, threads)

    def run_workflow_kma(self, fastq_input: FastqInput, sample_name: str, db_data: dict[str, Any],
                         threads: int = 8) -> None:
        """
        Runs the gene detection workflow using KMA.
        :param fastq_input: FASTQ input files
        :param sample_name: Sample name
        :param db_data: Database configuration
        :param threads: Number of threads to use
        :return: None
        """
        self.__create_input_fastq(fastq_input)
        config_data = self.__get_config_data(sample_name, db_data, 'kma', fastq_input.read_type)
        self.__run_workflow(config_data, threads)

    def __create_input_blast(self, fasta_path: Path) -> None:
        """
        Creates the input for the workflow in 'blast' mode.
        :param fasta_path: FASTA file path
        :return: None
        """
        path = Path(self._working_dir, gene_detection.INPUT_FASTA.format(db='db'))
        if not path.parent.exists():
            path.parent.mkdir(parents=True)
        snakemakeutils.dump_object([ToolIOFile(fasta_path)], path)

    def __create_input_fastq(self, fastq_input: FastqInput) -> None:
        """
        Creates the input for the workflow in FASTQ format.
        :param fastq_input: FASTQ input
        :return: None
        """
        path = self._working_dir / 'fq_dict.io'
        if fastq_input.is_pe:
            fq_dict = {'PE': fastq_input.pe}
        else:
            fq_dict = {'SE': fastq_input.se}
        snakemakeutils.dump_object(fq_dict, path)

    def __get_config_data(self, sample_name: str, db_data: dict[str, Any], detection_method: str, input_type: str)\
            -> dict[str, Any]:
        """
        Returns the configuration data for Snakemake.
        :param sample_name: Sample name
        :param db_data: Database information
        :param detection_method: Detection method
        :param input_type: Input type
        :return: Config data
        """
        return {
            'working_dir': str(self._working_dir),
            'sample_name': sample_name,
            'gene_detection': {
                'dbs': {'db': db_data},
                'options': {'method': detection_method}
            },
            'input_type': input_type
        }

    def __set_output(self, output_files: dict[str, str]) -> None:
        """
        Sets the output of the workflow.
        :param output_files: Output files dictionary
        :return: None
        """
        log_file_path = self._working_dir / 'camel.log'
        self._output = GeneDetectionOutput(
            report_section=snakemakeutils.load_object(Path(self._working_dir, output_files['report']))[0].value,
            detected_hits=[v.value for v in snakemakeutils.load_object(Path(self._working_dir, output_files['hits']))],
            informs=snakemakeutils.load_object(Path(self._working_dir, output_files['informs'])),
            log_file=log_file_path if log_file_path.exists() else None
        )

    @property
    def output(self) -> GeneDetectionOutput:
        """
        Returns the output generated by the gene detection workflow.
        :return: Output
        """
        return self._output
