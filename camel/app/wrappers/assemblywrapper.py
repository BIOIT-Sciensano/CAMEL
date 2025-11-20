from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Any, Union

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.scriptutils.basepipe.fastqinput import FastqInput
from camel.app.core.snakemake import snakepipelineutils, snakemakeutils
from camel.snakefiles import assembly


@dataclass
class AssemblyOutput:
    """
    Holds the assembly output data.
    """
    report_section: HtmlReportSection
    tsv_summary: Path
    fasta_contigs: Path
    log_file: Optional[Path]
    informs: list[dict[str, Any]]
    qc_stats: Optional[dict[str, Any]] = None


class AssemblyWrapper:
    """
    This class is used as a wrapper class around the assembly Snakemake workflow.
    """

    def __init__(self, working_dir: Path, input_type: str = 'illumina') -> None:
        """
        Initializes the read trimming helper.
        :param working_dir: Working directory
        :param input_type: Input type
        :return: None
        """
        self._working_dir = working_dir
        self._working_dir.mkdir(parents=True, exist_ok=True)
        self._output = None
        self._input_type = input_type

    @property
    def assembler_key(self) -> str:
        """
        Returns the key for the assembler.
        :return: Assembler key
        """
        if self._input_type == 'illumina':
            return 'spades'
        elif self._input_type == 'ont':
            return 'flye'
        raise ValueError(f'Invalid read type: {self._input_type}')

    def run(self, name: str, fastq_in: FastqInput, min_ctg_len: Optional[int] = None,
            assembler_opts: Optional[dict] = None, calc_qc_stats: bool = False, threads: int = 8) -> None:
        """
        Runs the assembly workflow for paired-end input.
        :param name: Dataset name
        :param fastq_in: FASTQ input
        :param min_ctg_len: Minimum contig length
        :param assembler_opts: Additional options for the assembler
        :param calc_qc_stats: If True, determines the QC stats
        :param threads: Number of threads
        :return: None
        """
        if fastq_in.is_pe:
            fq_dict = {'PE': fastq_in.pe}
            if fastq_in.se_fwd is not None:
                fq_dict['SE_FWD'] = fastq_in.se_fwd
            if fastq_in.se_rev is not None:
                fq_dict['SE_REV'] = fastq_in.se_rev
        else:
            fq_dict = {'SE': fastq_in.se}
        snakemakeutils.dump_object(fq_dict, self._working_dir / 'fq_dict.io')
        self.__run_workflow(name, min_ctg_len, assembler_opts, calc_qc_stats, threads)

    def __run_workflow(
            self, name: str, min_ctg_len: Optional[int] = None, assembler_opts: Optional[dict] = None,
            calc_qc_stats: bool = False, threads: int = 8) -> None:
        """
        Runs the underlying workflow.
        :return: None
        """
        # Creat working directory
        if not self._working_dir.exists():
            self._working_dir.mkdir(parents=True)

        # Create config file
        config_data = self.__get_config_data(name, min_ctg_len, assembler_opts)
        config_file = snakepipelineutils.generate_config_file(config_data, self._working_dir)

        # Collect output files
        output_files = self.__get_output_files_dict(config_data, min_ctg_len, calc_qc_stats)
        snakepipelineutils.run_snakemake(
            assembly.SNAKEFILE, config_file, list(output_files.values()), Path(self._working_dir), threads)
        self.__set_output(output_files)

    def __get_output_files_dict(
            self, config_data: dict[str, Any], min_ctg_len: Union[int, None], calc_qc_stats: bool) -> dict[str, Path]:
        """
        Returns the dictionary with output files.
        :param config_data: Configuration data
        :param min_ctg_len: Minimum contig length
        :param calc_qc_stats: If True, QC stats are calculated
        :return: Dictionary with output files.
        """
        output_files: dict[str, str | Path] = {
            'HTML': assembly.OUTPUT_REPORT,
            'TSV': assembly.OUTPUT_SUMMARY,
            'FASTA': assembly.OUTPUT_FASTA,
            **{f'INFORMS_assembler_{idx}': p for idx, p in enumerate(assembly.get_command_informs(config_data))}
        }
        if min_ctg_len is not None:
            output_files['INFORMS_seqtk'] = assembly.OUTPUT_INFORMS_FILTERING
        if calc_qc_stats is True:
            key_fq = 'fastq_pe' if self._input_type == 'illumina' else 'fastq_se'
            output_files['INFORMS_mapper'] = assembly.get_mapping_inform(key_fq)
            output_files['INFORMS_depth'] = assembly.get_depth_inform(key_fq)
        return output_files

    def __get_config_data(self, name: str, min_ctg_len: Union[int, None], assembler_opts: Optional[dict] = None) -> \
            dict[str, Any]:
        """
        Builds the configuration file to run the assembly workflow.
        :param name: Dataset name
        :param min_ctg_len: Minimal contig length
        :param assembler_opts: Additional options for the assembler
        :return: Config data
        """
        config_data: dict[str, Any] = {
            'working_dir': str(self._working_dir),
            'assembly': {self.assembler_key: {}},
            'input': {
                'type': self._input_type,
                'sample_name': name
            }
        }

        # Assembler options
        config_data['assembly'][self.assembler_key] = assembler_opts if assembler_opts is not None else {}

        # Length filtering
        if min_ctg_len is not None:
            config_data['assembly']['min_contig_length'] = min_ctg_len
        return config_data

    def __set_output(self, output_files: dict[str, Path]) -> None:
        """
        Runs the Snakemake workflow.
        :param output_files: Output files dictionary
        :return: None
        """
        log_file_path = self._working_dir / 'camel.log'
        informs = [
            snakemakeutils.load_object(self._working_dir / p)
            for key, p in output_files.items() if key.startswith('INFORMS_assembler_')
        ]
        if 'INFORMS_seqtk' in output_files:
            informs.append(snakemakeutils.load_object(self._working_dir / output_files['INFORMS_seqtk']))
        if all(key in output_files for key in ('INFORMS_mapper', 'INFORMS_depth')):
            qc_stats = {
                'depth': snakemakeutils.load_object(self._working_dir / output_files['INFORMS_depth']),
                'mapping': snakemakeutils.load_object(self._working_dir / output_files['INFORMS_mapper']),
            }
        else:
            qc_stats = None
        self._output = AssemblyOutput(
            report_section=snakemakeutils.load_object(self._working_dir/ output_files['HTML'])[0].value,
            tsv_summary=self._working_dir / output_files['TSV'],
            fasta_contigs=snakemakeutils.load_object(self._working_dir / output_files['FASTA'])[0].path,
            informs=informs,
            log_file=log_file_path if log_file_path.exists() else None,
            qc_stats=qc_stats
        )

    @property
    def output(self) -> AssemblyOutput:
        """
        Returns the output of the assembly workflow.
        :return: Assembly output
        """
        return self._output
