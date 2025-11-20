#!/usr/bin/env python
import dataclasses
import json
import shutil
from pathlib import Path
from typing import Any

import click

from camel.app.cli import cliutils
from camel.app.core.reports import reportutils
from camel.app.core.reports.htmlreport import HtmlReport
from camel.app.loggers import initialize_logging
from camel.app.scriptutils import inputhelper, model
from camel.app.scriptutils.basepipe import basepipeutils
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.app.scriptutils.inputhelper import helper_by_input_type, ONTHelper
from camel.app.scriptutils.inputhelper.inputhelperbase import InputHelperBase
from camel.app.scriptutils.model import BaseOptions
from camel.app.wrappers.genedetectionwrapper import (
    GeneDetectionWrapper,
    GeneDetectionOutput,
)


@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Custom options for the gene detection script.
    """
    db: Path = dataclasses.field(default=None, metadata={'help': 'Path to gene detection database'})
    working_dir: Path = dataclasses.field(default=Path.cwd(), metadata={'help': 'Working directory'})
    detection_method: str | None = dataclasses.field(default='blast', metadata={
        'help': 'Gene detection method',
        'choices': ['blast', 'kma'],
        'show_default': True})
    threads: int | None = dataclasses.field(default=1, metadata={'help': 'Number of threads', 'show_default': True})

    # BLAST options
    blast_min_percent_identity: int | None = dataclasses.field(default=90, metadata={
        'help': 'Minimum percent identity for BLAST',
        'show_default': True})
    blast_min_percent_coverage: int | None = dataclasses.field(default=60, metadata={
        'help': 'Minimum percent coverage for BLAST',
        'show_default': True})
    blast_task: str | None = dataclasses.field(default='megablast', metadata={
        'help': 'BLAST task', "show_default": True})
    blast_filtering_method: str = dataclasses.field(default='cluster', metadata={
        'help': 'BLAST filtering method',
        'choices': ['cluster', 'score', 'overlap'],
        'show_default': True
    })
    blast_score_nb_of_hits: int | None = dataclasses.field(default=5, metadata={
        'help': 'Number of hits to use for BLAST filtering (for score filtering only)',
        'show_default': True,
    })
    blast_reads: bool = dataclasses.field(default=False, metadata={
        'help': 'Perform BLAST search of the reads directly instead of on the assembly (ONT input)'})

    # KMA options
    kma_min_percent_identity: int | None = dataclasses.field(default=90, metadata={
        'help': 'Minimum percent identity for KMA',
        'show_default': True})
    kma_min_percent_coverage: int | None = dataclasses.field(default=60, metadata={
        'help': 'Minimum percent coverage for KMA',
        'show_default': True})
    kma_ont: bool = dataclasses.field(default=False, metadata={'help': 'Use ONT preset for KMA'})
    kma_cge: bool = dataclasses.field(default=False, metadata={'help': 'Use CGE preset for KMA'})
    kma_apm: str | None = dataclasses.field(default='u', metadata={
        'help': 'APM mode for KMA', 'choices': ['u', 'p', 'f']})


class MainGeneDetection(BaseScript[ScriptInput, ScriptOutput, Options]):
    """
    Main script for the gene detection tool.
    """

    def __init__(self, script_in: ScriptInput, helper: InputHelperBase, script_out: ScriptOutput, script_opts: Options) -> None:
        """
        Initializes the main script.
        :param script_in: Script input
        :param helper: Input helper
        :param script_out: Script output
        :param script_opts: Script options
        """
        super().__init__(
            name='Gene detection',
            version='1.0',
            script_in=script_in,
            script_out=script_out,
            script_opts=script_opts
        )
        self._helper: InputHelperBase = helper

    def _execute(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        # Initialize report
        report = reportutils.init_report(
            path_out=Path(self._script_out.html),
            key='Gene detection',
            title='Gene detection',
            dir_out=self._script_out.dir
        )
        report.add_html_object(reportutils.create_overview_section(
            version=self._version,
            dataset_name=self._script_in.name,
            input_file_str=self._script_in.input_str,
            extra_data=[('Detection method', self._script_opts.detection_method)]
        ))
        report.save()

        # Prepare wrapper
        wrapper = GeneDetectionWrapper(self._helper.working_dir)
        db_data = self.__get_db_metadata()

        # Run wrapper
        if self._script_opts.detection_method == 'blast':
            if self._script_opts.blast_reads:
                if not isinstance(self._helper, ONTHelper):
                    raise ValueError('BLAST reads is only supported for ONT input')
                fasta_input = self._helper.prepare_fastq_read_input(self._script_in, report)
            else:
                fasta_input = self._helper.prepare_fasta_input(self._script_in, report)
            # Save assembly if specified
            if self._script_out.fasta is not None:
                shutil.copyfile(str(fasta_input), self._script_out.fasta)
            wrapper.run_blast(fasta_input, self._script_in.name, db_data, self._script_opts.threads)
        elif self._script_opts.detection_method == 'kma':
            fastq_input = self._helper.prepare_fastq_input(self._script_in, report)
            wrapper.run_kma(fastq_input, self._script_in.name, db_data, self._script_opts.threads)

        # Export all output
        self.__export_output(report, wrapper.output)

    def __get_db_metadata(self) -> dict[str, Any]:
        """
        Returns the database information dictionary.
        :return: Database information dictionary
        """
        config_data = {'path': str(self._script_opts.db)}

        # Add specific options
        if self._script_opts.detection_method == 'blast':
            basepipeutils.dict_merge(
                config_data,
                {
                    'params': {
                        'blastn': {
                            'blast_reads': True if self._script_opts.blast_reads else False,
                            'filtering_method': self._script_opts.blast_filtering_method,
                            'min_coverage': self._script_opts.blast_min_percent_coverage,
                            'min_percent_identity': self._script_opts.blast_min_percent_identity,
                            'score_nb_of_hits': self._script_opts.blast_score_nb_of_hits,
                            'task': self._script_opts.blast_task
                        }
                    }
                }
            )
        elif self._script_opts.detection_method == 'kma':
            basepipeutils.dict_merge(
                config_data,
                {
                    'params': {
                        'kma': {
                            'min_percent_identity': self._script_opts.kma_min_percent_identity,
                            'min_coverage': self._script_opts.kma_min_percent_coverage,
                            'ont': self._script_opts.kma_ont,
                            'cge': self._script_opts.kma_cge,
                            'apm': self._script_opts.kma_apm
                        }
                    }
                }
            )

        # Add the extra metadata column
        with (self._script_opts.db / 'db_metadata.txt').open() as handle:
            db_metadata = json.load(handle)
            if 'extra_column' in db_metadata:
                config_data['metadata'] = db_metadata['extra_column']
        return config_data

    def __export_output(self, report: HtmlReport, output: GeneDetectionOutput) -> None:
        """
        Exports the output of the workflow.
        :param report: HTML report
        :param output: Workflow output
        :return: None
        """
        self._helper.logs['gene_detection'] = output.log_file if output.log_file is not None else None
        self._helper.informs.append(output.informs)
        self._helper.export_output_and_commands_section(report, output.report_section)


@click.command(name='gene_detection', short_help='Detection of genes in input FASTA / FASTQ files')
@basescriptutils.add_input_opts(supported=[model.InputType.FASTA, model.InputType.ILLUMINA, model.InputType.ONT])
@basescriptutils.add_output_opts
@cliutils.add_click_options_from_dataclass(Options)
@inputhelper.add_helper_opts
def main(**kwargs) -> None:
    """
    Gene detection.
    """
    # Parse the script input
    script_input = basescriptutils.parse_script_input(kwargs)
    script_opts = Options(**cliutils.from_kwargs(Options, kwargs))

    # Initialize the helper class to prepare the input
    helper = helper_by_input_type[script_input.type_](dir_=script_opts.working_dir, name=script_input.name)
    helper.set_opts(*helper.opts_from_cli(kwargs))

    # Run the main script
    script = MainGeneDetection(
        script_in=basescriptutils.parse_script_input(kwargs),
        script_out=basescriptutils.parse_script_output(kwargs),
        script_opts=script_opts,
        helper=helper
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
