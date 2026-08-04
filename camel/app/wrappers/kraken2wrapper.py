from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Union

from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.reports.htmlreportsection import HtmlReportSection

from camel.app.core.snakemake import snakemakeutils, snakepipelineutils
from camel.app.loggers import logger
from camel.app.scriptutils.basepipe.fastqinput import FastqInput
from camel.snakefiles import contamination_check_kraken


@dataclass
class Kraken2Output:
    """
    Dataclass to store the Kraken2 output.
    """
    report_section: HtmlReportSection
    tsv_summary: Path
    informs: dict[str, Any]
    informs_commands: list[dict[str, Any]]
    log_file: Optional[Path]


class Kraken2Wrapper:
    """
    This class is used as a wrapper class around the Kraken2 contamination check workflow.
    """

    def __init__(self, working_dir: Union[str, Path]) -> None:
        """
        Initializes the Kraken 2 wrapper.
        :param working_dir: Working directory
        :return: None
        """
        self._working_dir = Path(working_dir)
        self._output = None

    def run_fastq(
            self, sample_name: str, fastq_input: FastqInput, read_type: str, expected_species: str, db: Path,
            level_of_depth: str = 'S', threads: int = 8) -> None:
        """
        Runs the workflow on FASTQ input data.
        :param sample_name: Sample name
        :param fastq_input: FASTQ input
        :param read_type: Read type
        :param expected_species: Expected species
        :param db: Database
        :param level_of_depth: Species ('S') or Genus ('G') level of contamination check
        :param threads: Number of threads
        :return: None
        """
        snakemakeutils.dump_object(fastq_input.to_fq_dict(), self._working_dir / 'fq_dict.io')
        self.__run(read_type, sample_name, db, expected_species, level_of_depth, threads)

    def run_fasta(
            self, sample_name: str, fasta_in: Path, expected_species: str, db: Path, level_of_depth: str = 'S',
            threads: int = 8) -> None:
        """
        Runs the workflow on FASTA input data.
        :param sample_name: Sample name
        :param fasta_in: FASTA input
        :param expected_species: Expected species
        :param db: Database
        :param level_of_depth: Species ('S') or Genus ('G') level of contamination check
        :param threads: Number of threads
        :return: None
        """
        dir_fasta_in = self._working_dir / 'assembly' / 'filtering'
        dir_fasta_in.mkdir(parents=True, exist_ok=True)
        snakemakeutils.dump_object([ToolIOFile(fasta_in)], dir_fasta_in / 'fasta.io')
        self.__run('fasta', sample_name, db, expected_species, level_of_depth, threads)

    def __run(
            self, input_type: str, sample_name: str, db: Path, expected_species: str, level_of_depth: str = 'S',
            threads: int = 8) -> None:
        """
        Runs the kraken2 workflow.
        :param input_type: input type
        :param sample_name: Sample name
        :param expected_species: Expected species
        :param db: Database
        :param level_of_depth: Species ('S') or Genus ('G') level of contamination check
        :param threads: Number of threads
        :return: None
        """
        # Config file
        config_data = {
            'contamination_check': {
                'db': str(db),
                'expected_species': expected_species,
                'level_of_depth': level_of_depth},
            'working_dir': str(self._working_dir),
            'input': {'type': input_type, 'sample_name': sample_name}
        }
        config_file = snakepipelineutils.generate_config_file(config_data, self._working_dir)

        # Output files
        output_files: dict[str, str | Path] = {
            'HTML': contamination_check_kraken.OUTPUT_REPORT,
            'TSV': str(contamination_check_kraken.OUTPUT_SUMMARY).format(input_format='{input_format}', ext='tsv'),
            'INFORMS': contamination_check_kraken.OUTPUT_INFORMS,
            'INFORMS_KRAKEN': contamination_check_kraken.OUTPUT_INFORMS_KRAKEN2,
        }
        try:
            input_format = {'fasta': 'fasta', 'illumina': 'fastq_pe', 'ont': 'fastq_se'}[input_type]
        except KeyError as err:
            logger.error(f'Unsupported input type: {input_type}')
            raise err
        for k, p in output_files.items():
            if '{input_format}' not in str(p):
                continue
            output_files[k] = Path(str(p).format(input_format=input_format))

        # Run Snakemake
        snakepipelineutils.run_snakemake(
            contamination_check_kraken.SNAKEFILE, config_file, [Path(x) for x in output_files.values()],
            self._working_dir, threads)
        self.__set_output(output_files)

    def __set_output(self, output_files: dict[str, str]) -> None:
        """
        Runs the Snakemake workflow.
        :param output_files: Output files dictionary
        :return: None
        """
        log_file_path = self._working_dir / 'camel.log'
        self._output = Kraken2Output(
            report_section=snakemakeutils.load_object(self._working_dir / output_files['HTML'])[0].value,
            tsv_summary=self._working_dir / output_files['TSV'],
            informs=snakemakeutils.load_object(self._working_dir / output_files['INFORMS']),
            informs_commands=[snakemakeutils.load_object(self._working_dir / output_files['INFORMS_KRAKEN'])],
            log_file=log_file_path if log_file_path.exists() else None
        )

    @property
    def output(self) -> Kraken2Output:
        """
        Returns the output of the assembly workflow.
        :return: Assembly output
        """
        return self._output
