#!/usr/bin/env python
import dataclasses
from datetime import datetime
from pathlib import Path

import click
from camelcore.app.io.tooliodirectory import ToolIODirectory
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.utils import reportutils

from camel.app.cli import cliutils
from camel.app.loggers import initialize_logging
from camel.app.scriptutils import model
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.app.scriptutils.model import BaseOptions
from camel.app.tools.resfinder.resfinder import ResFinder
from camel.app.tools.resfinder.resfinderreporter import ResFinderReporter
from camel.resources import DIR_CITATIONS
from camel.version import __VERSION__


@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Specific options for ResFinder.
    """
    db: Path = dataclasses.field(metadata={'help': 'Path to ResFinder database'})
    working_dir: Path = dataclasses.field(default=Path.cwd(), metadata={'help': 'Working directory'})
    acquired: bool = dataclasses.field(default=False, metadata={'help': 'Use acquired genes'})
    point: bool = dataclasses.field(default=False, metadata={'help': 'Screen for point mutations'})
    min_cov: int = dataclasses.field(default=60, metadata={'help': 'Minimum coverage threshold'})
    threshold: int = dataclasses.field(default=80, metadata={'help': 'Minimum identity threshold'})
    acq_overlap: int = dataclasses.field(default=30, metadata={'help': 'Overlap for acquired genes'})
    species: str | None = dataclasses.field(default=None, metadata={'help': 'Species'})
    threads: int | None = dataclasses.field(default=1, metadata={'help': 'Number of threads', 'show_default': True})


class MainResFinder(BaseScript[ScriptInput, ScriptOutput, Options]):
    """
    This class is used to run the main ResFinder local script.
    """

    def __init__(self, in_: ScriptInput, out: ScriptOutput, opts: Options) -> None:
        """
        Initializes the main script.
        :param in_: Script input
        :param out: Script output
        :param opts: Options
        :return: None
        """
        tool_version = ResFinder().version
        super().__init__(
            name='ResFinder',
            version=f'{tool_version}+CAMEL_{__VERSION__}',
            script_in=in_,
            script_out=out,
            script_opts=opts
        )

    def _execute(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        # Initialize report
        report = reportutils.init_report(
            path_out=self._script_out.html,
            dir_out=self._script_out.dir,
            key=self.name,
            title=self.name,
        )
        additional_info = [
            ['Species', '<i>{}</i>'.format(self._script_opts.species.replace('"', '')) if
                self._script_opts.species is not None else 'Not specified'],
            ['Min % identity',
             f'{self._script_opts.threshold}' if self._script_opts.threshold is not None else 'Curated (default)'],
            ['Min % coverage', f'{self._script_opts.min_cov}'],
        ]
        report.add_html_object(reportutils.create_overview_section(
            version=self.version,
            dataset_name=self._script_in.name,
            input_file_str=self._script_in.input_str,
            date=datetime.now(),
            extra_data=additional_info
        ))
        report.save()

        # Run tools
        resfinder = self.__run_resfinder()
        section = self.__run_reporter(resfinder)
        report.add_html_object(section)
        section.copy_files(report.output_dir)

        # Save report
        all_informs = [resfinder.informs]
        report.add_html_object(reportutils.create_commands_section(all_informs, self._script_opts.working_dir))
        report.add_html_object(reportutils.create_citations_section(DIR_CITATIONS, ['Bortolaia_2020-resfinder_4.0']))
        report.save()

    def __run_resfinder(self) -> ResFinder:
        """
        Runs ResFinder.
        :return: ResFinder tool instance.
        """
        resfinder = ResFinder()
        # Input files
        resfinder.add_input_files({'DIR': [ToolIODirectory(self._script_opts.db)]})
        if self._script_in.type_ is model.InputType.FASTA:
            resfinder.add_input_files({'FASTA': [ToolIOFile(self._script_in.fasta)]})
        elif self._script_in.type_ is model.InputType.ILLUMINA:
            resfinder.add_input_files({'FASTQ_PE': [ToolIOFile(x) for x in self._script_in.fastq_pe]})
        elif self._script_in.fastq_se is not None:
            resfinder.add_input_files({'FASTQ_SE': [ToolIOFile(self._script_in.fastq_se)]})

        # Update parameters
        resfinder.update_parameters(
            min_cov=0.6,
            kma_threads=self._script_opts.threads,
            output_path=str(self._script_opts.working_dir),
            threshold=0.8)
        if self._script_opts.min_cov != 60:
            resfinder.update_parameters(min_cov=self._script_opts.min_cov / 100.0)
        if self._script_opts.threshold != 80:
            resfinder.update_parameters(threshold=self._script_opts.threshold / 100.0)
        if self._script_opts.point is True:
            try:
                resfinder.update_parameters(point=True, species='"' + self._script_opts.species.replace('_', ' ') + '"')
            except AttributeError:
                raise ValueError('--point requires a --species argument')
        if self._script_opts.acquired is not None:
            resfinder.update_parameters(acquired=True, acq_overlap=self._script_opts.acq_overlap)

        # Run the tools
        resfinder.run(self._script_opts.working_dir)
        return resfinder

    def __run_reporter(self, resfinder: ResFinder) -> HtmlReportSection:
        """
        Runs resfinder reporter.
        :param resfinder: ResFinder tool instance.
        :return: None.
        """
        reporter = ResFinderReporter()
        reporter.add_input_files({'TSV_pheno_general': resfinder.tool_outputs['TSV_pheno_general']})
        if self._script_opts.acquired is True:
            reporter.add_input_files({'TSV_genes': resfinder.tool_outputs['TSV_genes']})
        if self._script_opts.point is True:
            reporter.add_input_files({
                'TSV_point': resfinder.tool_outputs['TSV_point'],
                'TSV_pheno_species': resfinder.tool_outputs['TSV_pheno_species']})
        reporter.add_input_informs({'resfinder': resfinder.informs})
        reporter.run()
        return reporter.tool_outputs['VAL_HTML'][0].value


@click.command(name='resfinder', short_help='Wrapper for ResFinder4')
@basescriptutils.add_input_opts()
@basescriptutils.add_output_opts
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Wrapper for ResFinder4.
    """
    script = MainResFinder(
        in_=basescriptutils.parse_script_input(kwargs),
        out=basescriptutils.parse_script_output(kwargs),
        opts=Options(**cliutils.from_kwargs(Options, kwargs))
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
