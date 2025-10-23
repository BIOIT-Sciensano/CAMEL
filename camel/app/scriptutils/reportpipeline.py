import abc
import argparse
import shutil
from pathlib import Path
from typing import Any, Union

from camel.app.core.reports import reportutils
from camel.app.scriptutils.basepipeline import BasePipeline
from camel.app.scriptutils.fastqinput import FastqInput
from camel.app.core.utils import fastautils
from camel.app.scriptutils import mainscriptutils, absolute_path_by_pathlib
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import assembly


class ReportPipeline(BasePipeline, metaclass=abc.ABCMeta):
    """
    Baseclass for pipelines with a report (HTML) and tabular (TSV) output.
    """

    @staticmethod
    def add_common_arguments(argument_parser: argparse.ArgumentParser) -> None:
        """
        Adds the common arguments for the pipeline.
        :param argument_parser: Argument parser
        :return: None
        """
        BasePipeline.add_common_arguments(argument_parser)

        # Output
        argument_parser.add_argument(
            '--output-dir',
            help='Output directory',
            type=absolute_path_by_pathlib,
            default=Path(Path.cwd(), 'out'))
        argument_parser.add_argument(
            '--output-html',
            type=absolute_path_by_pathlib,
            default=Path(Path.cwd(), 'out', 'report.html'))
        argument_parser.add_argument(
            '--output-tsv',
            help='Output file for the summary in TSV format',
            type=absolute_path_by_pathlib,
            default=Path(Path.cwd(), 'out', 'summary.tsv'))
        argument_parser.add_argument(
            '--output-json',
            help='Output file for the summary JSON format',
            type=absolute_path_by_pathlib)
        argument_parser.add_argument(
            '--output-fasta',
            type=absolute_path_by_pathlib,
            help='output path for assembled contigs')

        # Options
        argument_parser.add_argument(
            '--detection-method', help="Method for allele / gene detection", choices=['blast', 'kma', 'mist'], default='blast')
        argument_parser.add_argument(
            '--report-include-fastq', help="Include the FASTQ files in the report", action='store_true')
        argument_parser.add_argument(
            '--report-include-bam', help="Include the BAM file in the report", action='store_true')

        # Parameters
        argument_parser.add_argument(
            '--cov-max', default=100.0, type=float,
            help='Maximum coverage (datasets with higher estimated coverage will be downsampled to the given value)')

    def get_template_data(self, dict_input: dict) -> dict[str, Any]:
        """
        Returns the template data that is common to all pipelines.
        :param dict_input: Dictionary with pipeline input files
        :return: Template data
        """
        # Add common entries
        config_data = super().get_template_data(dict_input)

        # Add report-specific entries
        mainscriptutils.dict_merge(config_data, {
            'output_dir': str(self._args.output_dir),
            'output_report': str(self._args.output_html),
            'output_tabular': str(self._args.output_tsv),
            'output_json': str(self._args.output_json) if self._args.output_json is not None else None,
            'sequence_typing': {'options': {'method': self._args.detection_method}},
            'gene_detection': {'options': {'method': self._args.detection_method}},
        })

        # FASTQ export
        if self._args.report_include_fastq is True:
            config_data['read_trimming']['export_fastq'] = True

        # Technology-specific options
        if (self._args.input_type == 'illumina') and (self._args.library is not None):
            config_data['read_trimming']['adapter'] = self._args.library

        # Empty variant filtering config
        config_data['variant_filtering'] = {}

        # Human read scrubbing
        if 'human_read_scrubbing' in self._args:
            config_data['read_scrubbing'] = {}

        return config_data

    def _validate_input_files(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        mainscriptutils.validate_input_files(self._args)

    def _export_assembly(self) -> None:
        """
        Exports the assembly to the specified output location (optional).
        :return: None
        """
        if self._args.output_fasta is None:
            logger.debug('Not exporting assembly')
            return
        path_io = self._args.working_dir / assembly.OUTPUT_FASTA
        path_fasta = snakemakeutils.load_object(path_io)[0].path
        fastautils.rename_sequences_regex(path_fasta, self._args.output_fasta, '', '', description=self.sample_name)
        logger.info(f'Output FASTA file exported to: {self._args.output_fasta}')

    @staticmethod
    def format_input_string(dict_in: dict[str, Any]) -> str:
        """
        Formats the input string based on the dictionary with pipeline input files.
        :param dict_in: Pipeline input dictionary
        :return: Formatted input string
        """
        return ', '.join([entry['name'] for key, entries in dict_in.items() for entry in entries])

    @staticmethod
    def construct_fq_dict(snake_in: Any, input_type: str, path_out: Path) -> None:
        """
        Constructs a dictionary with the FASTQ input based on the snakemake input.
        :param snake_in: Snakemake input
        :param input_type: Input type
        :param path_out: Path to store dictionary
        :return: None
        """
        # FASTA input
        if input_type in ('fasta', 'fasta_with_vcf'):
            fq_pe = snakemakeutils.load_object(Path(snake_in.FASTQ_PE))
            snakemakeutils.dump_object(FastqInput('illumina', fq_pe, is_pe=True).to_fq_dict(), path_out)

        # PE reads (illumina)
        elif input_type == 'illumina':
            shutil.copyfile(snake_in.FASTQ_PE, path_out)

        # SE reads (iontorrent, ont)
        elif input_type in ('iontorrent', 'ont'):
            shutil.copyfile(snake_in.FASTQ_SE, path_out)

        # Hybrid reads
        elif input_type == 'hybrid':
            fq_pe = FastqInput.from_fq_dict(Path(snake_in.FASTQ_PE), 'illumina')
            fq_se = FastqInput.from_fq_dict(Path(snake_in.FASTQ_SE), 'ont')
            snakemakeutils.dump_object(FastqInput(
                'hybrid', pe=fq_pe.pe, se_fwd=fq_pe.se_fwd, se_rev=fq_pe.se_rev, se=fq_se.se, is_pe=True,
                is_trimmed=True).to_fq_dict(), path_out)
        else:
            raise ValueError(f'Invalid input type: {input_type}')
        logger.info(f'FASTQ dict object created: {path_out}')

    @staticmethod
    def export_command_section(snake_in: Any, path_out: Path, dir_: Path) -> None:
        """
        Generates the command report section and exports it to the target path.
        :param snake_in: Snakemake input
        :param path_out: Output path
        :param dir_: Working directory (used to mask absolute paths)
        :return :None
        """
        informs = []
        for content in [snakemakeutils.load_object(Path(io)) for io in snake_in]:
            if type(content) is dict:
                informs.append(content)
            elif type(content) is list:
                informs.extend(content)
        section = reportutils.create_commands_section(informs, dir_)
        snakemakeutils.dump_object([ToolIOValue(section)], path_out)
        logger.info(f'Command exported to: {path_out}')

    @staticmethod
    def add_content_scrubbing(
            structure: list[tuple], input_type: str, reports_scrubbing: list[Union[Path, str]]) -> None:
        """
        Adds the report content for the human read scrubbing.
        :param structure: Report structure
        :param input_type: Input type
        :param reports_scrubbing: Human read scrubbing output report(s)
        :return: None
        """
        # Create dictionaries with the technology as key and the reports as values
        report_scrubbing_by_input_format = {
            p_html.parents[1].name: p_html for p_html in [Path(x) for x in reports_scrubbing]}

        # Add the report content
        if input_type in ('fasta', 'fasta_with_vcf'):
            structure.append(
                ('Human read removal', 'human read removal', [
                    report_scrubbing_by_input_format['fasta']]))
        elif input_type == 'illumina':
            structure.append(
                ('Human read removal', 'human read removal', [
                    report_scrubbing_by_input_format['fastq_pe']]))
        elif input_type == 'ont':
            structure.append(
                ('Human read removal', 'human read removal', [report_scrubbing_by_input_format['fastq_se']]))
        elif input_type == 'hybrid':
            structure.append(
                ('Human read removal - Illumina', 'human_read_removal_ilmn',
                 [report_scrubbing_by_input_format['fastq_pe']]))
            structure.append(
                ('Human read removal - ONT', 'human_read_removal_ont',
                 [report_scrubbing_by_input_format['fastq_se']]))

        else:
            raise ValueError(f'Invalid input type: {input_type}')

    @staticmethod
    def add_content_trim_basic_qc(
            structure: list[tuple], input_type: str, reports_ds: list[Union[Path, str]],
            reports_trim: list[Union[Path, str]]) -> None:
        """
        Adds the report content for the downsampling, basic QC, and read trimming.
        :param structure: Report structure
        :param input_type: Input type
        :param reports_ds: Downsampling output reports
        :param reports_trim: Trimming output reports
        :return: None
        """
        # Create dictionaries with the technology as key and the reports as values
        report_trim_by_tech = {
            next(p.name.replace('trimming_', '') for p in p_html.parents if p.name.startswith('trimming_')):
                p_html for p_html in [Path(x) for x in reports_trim]}
        report_ds_by_read_key = {
            p_html.parents[1].name: p_html for p_html in [Path(x) for x in reports_ds]}

        # Add the report content
        if input_type in ('fasta', 'fasta_with_vcf'):
            pass
        elif input_type == 'illumina':
            structure.append(('Read trimming and basic QC', 'trim', [
                report_ds_by_read_key['fastq_pe'], report_trim_by_tech['illumina']]))
        elif input_type == 'ont':
            structure.append(('Read trimming and basic QC', 'trim', [
                report_ds_by_read_key['fastq_se'], report_trim_by_tech['ont']]))
        elif input_type == 'hybrid':
            structure.append(('Read trimming and basic QC - Illumina', 'trim_ilmn', [
                report_ds_by_read_key['fastq_pe'], report_trim_by_tech['illumina']]))
            structure.append(('Read trimming and basic QC - ONT', 'trim_ont', [
                report_ds_by_read_key['fastq_se'], report_trim_by_tech['ont']]))

    @staticmethod
    def add_content_contamination_check(
            structure: list[tuple], input_type: str, reports_contamination: list[Union[Path, str]],
            report_confindr: Union[Path, str, None]) -> None:
        """
        Adds the report content for the contamination check.
        :param structure: Report structure
        :param input_type: Input type
        :param reports_contamination: Contamination check output report(s)
        :param report_confindr: ConFindr report
        :return: None
        """
        # Create dictionaries with the technology as key and the reports as values
        report_k2_by_input_format = {
            p_html.parents[1].name: p_html for p_html in [Path(x) for x in reports_contamination]}

        if input_type in ('fasta', 'fasta_with_vcf'):
            # FASTA input -> only Kraken2
            structure.append(
                ('Contamination check', 'contamination', [report_k2_by_input_format['fasta']]))
        elif input_type == 'illumina':
            reports = [report_k2_by_input_format['fastq_pe']]
            if report_confindr is not None:
                reports.append(Path(report_confindr))
            structure.append(('Contamination check', 'contamination', reports))
        elif input_type == 'ont':
            # ONT input -> Kraken2 and ConFindr
            reports = [report_k2_by_input_format['fastq_se']]
            if report_confindr is not None:
                reports.append(Path(report_confindr))
            structure.append(('Contamination check', 'contamination', reports))
        elif input_type == 'hybrid':
            structure.append(
                ('Contamination check', 'contamination',
                 [report_k2_by_input_format['fastq_pe'], report_k2_by_input_format['fastq_se'], Path(report_confindr)]))
        else:
            raise ValueError(f'Invalid input type: {input_type}')
