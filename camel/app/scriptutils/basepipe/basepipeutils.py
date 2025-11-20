import datetime
import gzip
import shutil
import socket
from pathlib import Path
from collections.abc import Mapping
from typing import Any, Union

from camel.app.config import config
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.reports import reportutils
from camel.app.core.snakemake import snakemakeutils
from camel.app.core.utils import fileutils
from camel.app.loggers import logger
from camel.app.scriptutils.basepipe.fastqinput import FastqInput

TIMESTAMP_FILENAME = "%Y%m%d-%H%M%S"


def get_timestamp_str(timestamp: datetime.datetime = datetime.datetime.now()) -> str:
    """
    Returns the given time stamp as a string that can be used in a filename.
    :param timestamp: Timestamp (default to current time)
    :return: Timestamp as string
    """
    return timestamp.strftime(TIMESTAMP_FILENAME)


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

def dict_merge(dct: dict[str, Any], merge_dct: dict) -> None:
    """
    Recursive dict merge. Inspired by :meth:``dict.update()``, instead of updating only top-level keys,
    dict_merge recurses down into dicts nested to an arbitrary depth, updating keys. The ``merge_dct`` is merged into
    ``dct`` (https://gist.github.com/angstwad/bf22d1822c38a92ec0a9).
    :param dct: dict onto which the merge is executed
    :param merge_dct: dct merged into dct
    :return: None
    """
    for k, v in merge_dct.items():
        if k in dct and isinstance(dct[k], dict) and isinstance(merge_dct[k], Mapping):
            dict_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]

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

def store_log_file(log_file: Path, basename: str, is_error_log: bool = False, dir_: Path | None = None) -> Path | None:
    """
    Stores a log file on disk.
    :param log_file: Log file to store
    :param basename: Basename
    :param is_error_log: Boolean to indicate if this is an error log
    :param dir_: (Optional) Directory to store file, defaults to value from config
    :return: Path to log file
    """
    if not is_error_log:
        raise ValueError("Only error logs can be stored.")

    dir_out = dir_ if dir_ is not None else config.dir_logs
    if dir_out is None:
        logger.warning('No directory specified for storing error logs, not copying log file')
        return None
    if not dir_out.exists():
        raise RuntimeError(f'Logging directory does not exist: {dir_out}')

    # Determine output file
    prefix = 'error' if is_error_log else 'camel'
    output_path = dir_ / '{}.txt.gz'.format('__'.join([
        prefix,
        fileutils.make_valid(basename).lower(),
        fileutils.make_valid(socket.gethostname()),
        get_timestamp_str()
    ]))

    # Copy the log file
    with gzip.open(output_path, 'wb') as file_out, log_file.open('rb') as file_in:
        # noinspection PyTypeChecker
        shutil.copyfileobj(file_in, file_out)
    logger.debug(f"Log file stored in: '{output_path}'")
    return output_path

def prepare_galaxy_output(output_dir: Path, output_html: Path) -> None:
    """
    Prepares the Galaxy output files at the start of the script.
    - The output HTML file is removed, so Snakemake can regenerate it
    - The output directory is created if it does not exist yet.
    :param output_dir: Output directory
    :param output_html: Output report path
    :return: None
    """
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    if output_html.exists():
        output_html.unlink()
