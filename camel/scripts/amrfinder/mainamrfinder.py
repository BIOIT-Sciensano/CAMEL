#!/usr/bin/env python
import dataclasses
import shutil
from datetime import datetime
from pathlib import Path

import click
from camelcore.app.io.tooliodirectory import ToolIODirectory
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.reports.htmlreport import HtmlReport
from camelcore.app.utils import fileutils, reportutils

from camel.app.cli import cliutils
from camel.app.loggers import initialize_logging
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.basescript.fastainput import FastaInput
from camel.app.scriptutils.model import BaseOptions, BaseOutput
from camel.app.tools.amrfinder.amrfinder import AMRFinder
from camel.app.tools.amrfinder.amrfinderreporter import AMRFinderReporter
from camel.resources import DIR_CITATIONS


@dataclasses.dataclass(frozen=True)
class Output(BaseOutput):
    """
    Tool output.
    """
    output_html: Path = dataclasses.field(metadata={'help': 'Output HTML file'})
    output_dir: Path | None = dataclasses.field(metadata={'help': 'Output directory'})
    output_tsv: Path | None = dataclasses.field(default=None, metadata={'help': 'Output TSV file'})

@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Specific options for AMFinder+.
    """
    db: Path = dataclasses.field(metadata={'help': 'Path to AMFinder+ database'})
    working_dir: Path = dataclasses.field(default=Path('working'), metadata={
        'help': 'Working directory', 'show_default': True})
    min_id: int | None = dataclasses.field(default=None, metadata={
        'help': 'Minimum identity threshold (defaults to curated cutoffs if not specified)'})
    min_cov: int | None = dataclasses.field(default=50, metadata={
        'help': 'Minimum coverage threshold', 'show_default': True})
    organism: str | None = dataclasses.field(default=None, metadata={'help': 'Organism'})
    threads: int = dataclasses.field(default=1, metadata={'help': 'Number of threads', 'show_default': True})


class MainAMRFinder(BaseScript[FastaInput, Output, Options]):
    """
    Main script for the AMRFinder tool.
    """

    def __init__(self, in_: FastaInput, out: Output, opts: Options) -> None:
        """
        Initializes the main script.
        :param in_: Script input
        :param out: Script output
        :param opts: Options
        :return: None
        """
        super().__init__(
            name='AMRFinder+',
            version='1.0.0',
            title='AMRFinder+',
            script_in=in_,
            script_out=out,
            script_opts=opts
        )

    def __initialize_report(self) -> HtmlReport:
        """
        Initializes the HTML output report.
        :return: HTML report
        """
        # Initialize report
        report = reportutils.init_report(
            path_out=self._script_out.output_html,
            dir_out=self._script_out.output_dir,
            key=self.name,
            title=f'{self.name} v{self.version}',
        )

        # Create the overview section
        additional_info = [
            ['Organism', self._script_opts.organism if self._script_opts.organism is not None else 'Not specified'],
            ['Min % identity', self._script_opts.min_id if self._script_opts.min_id is not None else 'Curated (default)'],
            ['Min % coverage', self._script_opts.min_cov],
        ]
        report.add_html_object(reportutils.create_overview_section(
            version=self.version,
            dataset_name=self._script_in.name,
            input_file_str=self._script_in.input_str,
            date=datetime.now(),
            extra_data=additional_info
        ))
        report.save()
        return report

    def _execute(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        report = self.__initialize_report()

        # Run AMRFinder
        amrfinder = AMRFinder()
        amrfinder.add_input_files({
            'FASTA': [ToolIOFile(self._script_in.fasta)],
            'DIR': [ToolIODirectory(self._script_opts.db)]
        })
        amrfinder.update_parameters(
            min_cov=self._script_opts.min_cov / 100.0,
            output_path=f'amrfinder_{fileutils.make_valid(self._script_in.name)}.tsv',
            threads=self._script_opts.threads
        )
        if self._script_opts.min_id is not None:
            amrfinder.update_parameters(min_ident=self._script_opts.min_id / 100.0)
        if self._script_opts.organism is not None:
            amrfinder.update_parameters(organism=self._script_opts.organism)
        self._script_opts.working_dir.mkdir(parents=True, exist_ok=True)
        amrfinder.run(self._script_opts.working_dir.absolute())

        # Create the output section
        amrfinder_reporter = AMRFinderReporter()
        amrfinder_reporter.add_input_files({'TSV': amrfinder.tool_outputs['TSV']})
        amrfinder_reporter.add_input_informs({'amrfinder': amrfinder.informs})
        amrfinder_reporter.run(self._script_opts.working_dir.absolute())
        report.add_html_object(amrfinder_reporter.tool_outputs['HTML'][0].value)
        amrfinder_reporter.tool_outputs['HTML'][0].value.copy_files(report.output_dir)

        # Add commands and citations
        report.add_html_object(reportutils.create_commands_section([amrfinder.informs], self._script_opts.working_dir.absolute()))
        report.add_html_object(
            reportutils.create_citations_section(
                dir_=DIR_CITATIONS,
                keys_other=['Feldgarden_2019-ndaro'],
        ))
        report.save()

        # Copy the TSV output file if specified
        if self._script_out.output_tsv is not None:
            shutil.copyfile(amrfinder.tool_outputs['TSV'][0].path, self._script_out.output_tsv)


@click.command(name='amrfinder', short_help='Wrapper for NCBI AMRFinder+')
@cliutils.add_click_options_from_dataclass(FastaInput)
@cliutils.add_click_options_from_dataclass(Output)
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Wrapper for NCBI AMRFinder+.
    """
    script = MainAMRFinder(
        in_=FastaInput(**cliutils.from_kwargs(FastaInput, kwargs)),
        out=Output(**cliutils.from_kwargs(Output, kwargs)),
        opts=Options(**cliutils.from_kwargs(Options, kwargs))
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
