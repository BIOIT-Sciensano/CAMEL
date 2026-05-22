#!/usr/bin/env python
import dataclasses
import shutil
from datetime import datetime
from pathlib import Path

import click
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import fileutils, reportutils

from camel.app.cli import cliutils
from camel.app.loggers import initialize_logging
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.basescript.fastainput import FastaInput
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.app.scriptutils.model import BaseOptions
from camel.app.tools.integronfinder.integronfinder import IntegronFinder
from camel.app.tools.integronfinder.integronfinderreporter import IntegronFinderReporter
from camel.resources import DIR_CITATIONS


@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Specific options for IntegronFinder.
    """
    working_dir: Path = dataclasses.field(default=Path('working'), metadata={
        'help': 'Working directory', 'show_default': True})
    local_max: bool = dataclasses.field(default=False, metadata={
        'help': 'Allows thorough local detection (slower but more sensitive).'})
    threads: int = dataclasses.field(default=4, metadata={'help': 'Number of threads to use'})


class MainIntegronFinder(BaseScript[FastaInput, ScriptOutput, Options]):
    """
    Main script for the IntegronFinder tool.
    """

    def __init__(self, in_: FastaInput, out: ScriptOutput, opts: Options) -> None:
        """
        Initializes the main script.
        :param in_: Script input
        :param out: Script output
        :param opts: Options
        :return: None
        """
        super().__init__(
            name='IntegronFinder',
            version='1.0.0',
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
        report.add_html_object(reportutils.create_overview_section(
            version=self.version,
            dataset_name=self._script_in.name,
            input_file_str=self._script_in.input_str,
            date=datetime.now(),
        ))
        report.save()

        # Run IntegronFinder
        integron_finder = IntegronFinder()
        integron_finder.add_input_files({'FASTA': [ToolIOFile(self._script_in.fasta)]})
        integron_finder.update_parameters(
            local_max=self._script_opts.local_max, threads=self._script_opts.threads)
        integron_finder.run(self._script_opts.working_dir)

        # Create the output section
        reporter = IntegronFinderReporter()
        reporter.add_input_files({'TSV': integron_finder.tool_outputs['TSV']})
        reporter.add_input_informs({'integron_finder': integron_finder.informs})
        reporter.update_parameters(name=fileutils.make_valid(self._script_in.name))
        reporter.run(self._script_opts.working_dir)
        report.add_html_object(reporter.tool_outputs['HTML'][0].value)
        reporter.tool_outputs['HTML'][0].value.copy_files(report.output_dir)

        # Add commands and citations
        report.add_html_object(reportutils.create_commands_section(
            [integron_finder.informs], self._script_opts.working_dir))
        report.add_html_object(reportutils.create_citations_section(
            dir_=DIR_CITATIONS, keys_other=['Neron_2022-integronfinder']))
        report.save()

        # Copy the TSV output file when specified
        if self._script_out.tsv is not None:
            shutil.copyfile(integron_finder.tool_outputs['TSV'][0].path, self._script_out.tsv)


@click.command(name='integronfinder', short_help='Wrapper for IntegronFinder')
@cliutils.add_click_options_from_dataclass(FastaInput)
@basescriptutils.add_output_opts
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Wrapper for IntegronFinder.
    """
    script = MainIntegronFinder(
        in_=FastaInput(**cliutils.from_kwargs(FastaInput, kwargs)),
        out=basescriptutils.parse_script_output(kwargs),
        opts=Options(**cliutils.from_kwargs(Options, kwargs))
    )
    script.run()



if __name__ == '__main__':
    initialize_logging()
    main()
