import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.reports.htmlreportsection import HtmlReportSection

from camel.app.core.snakemake import snakemakeutils, snakepipelineutils
from camel.app.loggers import logger
from camel.app.scriptutils.basepipe.fastqinput import FastqInput
from camel.app.toolkits.sequencetyping.sequencetypingutils import SequenceTypingUtils
from camel.snakefiles import sequence_typing


@dataclass
class SequenceTypingInput:
    """
    This class is used to construct the input for the sequence typing workflow.
    """
    sample_name: str
    db_path: Path
    input_type: str
    fasta: Optional[ToolIOFile] = None
    fastq: Optional[FastqInput] = None
    db_key: str = 'mlst'


@dataclass
class SequenceTypingOutput:
    """
    This class contains the output of the sequence typing workflow.
    """
    report_section: HtmlReportSection
    tsv: Optional[Path]
    informs: list[dict[str, Any]]
    log_file: Optional[Path] = None


class SequenceTypingWrapper:
    """
    This class is used as a wrapper class around the sequence typing Snakemake workflow.
    """

    def __init__(self, working_dir: Path) -> None:
        """
        Initializes the read sequence typing helper.
        :param working_dir: Working directory
        """
        self._working_dir = working_dir
        self._output = None

    def run_workflow_blast(
            self, workflow_input: SequenceTypingInput, blast_task: Optional[str] = None,
            threads: Optional[int] = 8) -> None:
        """
        Runs the BLAST-based sequence typing workflow.
        :param workflow_input: Workflow input
        :param blast_task: blast task
        :param threads: Number of threads to use
        :return: None
        """
        logger.info("Running sequence typing workflow (blast)")
        if not self._working_dir.exists():
            self._working_dir.mkdir(parents=True)
        if workflow_input.fasta is None:
            raise ValueError("FASTA input is required to run the blast analysis")
        self.__create_fasta_input(workflow_input.fasta.path, workflow_input.db_key)
        options = {'blastn_task': blast_task} if blast_task is not None else None
        config_path = self.__create_config_file(
            sample_name=workflow_input.sample_name,
            detection_method='blast',
            input_type=workflow_input.input_type,
            db_key=workflow_input.db_key,
            db_path=workflow_input.db_path,
            additional_options=options
        )
        self.__run_snakefile(config_path, workflow_input, threads)

    def run_workflow_mist(self, workflow_input: SequenceTypingInput, threads: Optional[int] = 8) -> None:
        """
        Runs the MiST-based sequence typing workflow.
        :param workflow_input: Workflow input
        :param threads: Number of threads to use
        :return: None
        """
        if not self._working_dir.exists():
            self._working_dir.mkdir(parents=True)
        if workflow_input.fasta is None:
            raise ValueError("FASTA input is required to run the blast analysis")
        self.__create_fasta_input(workflow_input.fasta.path, workflow_input.db_key)
        config_path = self.__create_config_file(
            sample_name=workflow_input.sample_name,
            detection_method='mist',
            input_type=workflow_input.input_type,
            db_key=workflow_input.db_key,
            db_path=workflow_input.db_path,
            additional_options={},
        )
        self.__run_snakefile(config_path, workflow_input, threads)

    def run_workflow_kma(self, workflow_input: SequenceTypingInput, threads: Optional[int] = 8) -> None:
        """
        Runs the KMA-based sequence typing workflow.
        :param workflow_input: Input for the workflow
        :param threads: Number of threads to use
        :return: None
        """
        if not self._working_dir.exists():
            self._working_dir.mkdir()
        if workflow_input.fastq is None:
            raise ValueError('FASTQ input is required to run the KMA analysis')
        scheme_metadata = SequenceTypingUtils.parse_scheme_metadata(workflow_input.db_path)
        if len([locus for locus in scheme_metadata['loci'] if locus['type'] == 'peptide']) > 0:
            if workflow_input.fasta is None:
                raise ValueError('FASTA input is required when there are protein loci in the scheme')
            else:
                self.__create_fasta_input(workflow_input.fasta.path, workflow_input.db_key)
        snakemakeutils.dump_object(workflow_input.fastq.to_fq_dict(), self._working_dir / 'fq_dict.io')
        config_path=self.__create_config_file(
            sample_name=workflow_input.sample_name,
            detection_method='kma',
            input_type=workflow_input.input_type,
            db_key=workflow_input.db_key,
            db_path=workflow_input.db_path,
            additional_options={}
        )
        self.__run_snakefile(config_path, workflow_input, threads)

    def __create_config_file(
            self, sample_name: str, detection_method: str, input_type: str, db_key: str, db_path: Path,
            additional_options: dict[str, Any]) -> str:
        """
        Creates the configuration file for the sequence typing workflow.
        :param sample_name: Sample name
        :param detection_method: Detection method
        :param input_type: Input type
        :param db_key: Database key (e.g. 'mlst')
        :param db_path: Database path
        :param additional_options: Additional options to add to the configuration file.
        :return: Config file path
        """
        data = {
            'working_dir': str(self._working_dir),
            'input': {'type': input_type, 'sample_name': sample_name},
            'sequence_typing': {
                'dbs': {
                    db_key: {
                        'path': str(db_path),
                        **additional_options
                    }},
                'options': {
                    'method': detection_method
                }
            }
        }
        return snakepipelineutils.generate_config_file(data, self._working_dir)

    def __create_fasta_input(self, fasta_path: Path, scheme_key: str) -> None:
        """
        Creates the FASTA input for the workflow.
        :param scheme_key: Scheme key
        :param fasta_path: FASTA file path
        :return: None
        """
        fasta_in = self._working_dir / str(sequence_typing.INPUT_FASTA).format(scheme=scheme_key)
        if not fasta_in.parent.exists():
            fasta_in.parent.mkdir(parents=True)
        logger.debug(f'Creating FASTA input: {fasta_in}')
        snakemakeutils.dump_object([ToolIOFile(fasta_path)], fasta_in)

    def __run_snakefile(self, config_path: str, workflow_input: SequenceTypingInput, threads: int) -> None:
        """
        Runs the Snakefile.
        :param config_path: Path to the configuration file
        :param workflow_input: Workflow input
        :param threads: Number of threads to use
        :return: None
        """
        # Create dictionary with output files
        output_files = {
            'report': sequence_typing.OUTPUT_REPORT.format(scheme=workflow_input.db_key),
            'informs': sequence_typing.OUTPUT_INFORMS.format(scheme=workflow_input.db_key),
        }

        # Check if all loci are of the same type (DNA, peptide), otherwise tabular output cannot be created
        with (workflow_input.db_path / 'scheme_metadata.json').open() as handle:
            data_scheme = json.load(handle)
        locus_type = data_scheme['loci'][0]['type']
        if all(locus['type'] == locus_type for locus in data_scheme['loci']):
            output_files['tsv'] = sequence_typing.OUTPUT_TSV.format(scheme=workflow_input.db_key, locus_type=locus_type)

        # Run Snakemake
        snakepipelineutils.run_snakemake(
            sequence_typing.SNAKEFILE, config_path, [Path(x) for x in output_files.values()], self._working_dir,
            threads)

        # Collect output
        log_file_path = self._working_dir / 'camel.log'
        self._output = SequenceTypingOutput(
            report_section=snakemakeutils.load_object(self._working_dir / output_files['report'])[0].value,
            tsv=snakemakeutils.load_object(self._working_dir / output_files['tsv'])[0].path if 'tsv' in output_files else None,
            informs=snakemakeutils.load_object(self._working_dir / output_files['informs']),
            log_file=log_file_path if log_file_path.exists() else None
        )

    @property
    def output(self) -> SequenceTypingOutput:
        """
        Returns the output generated by the gene detection workflow.
        :return: Output
        """
        return self._output
