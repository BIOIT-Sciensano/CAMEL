#!/usr/bin/env python
import dataclasses
from pathlib import Path

import click
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import reportutils

from camel.app.cli import cliutils
from camel.app.loggers import initialize_logging
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.basescript.fastainput import FastaInput
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.app.scriptutils.model import BaseOptions
from camel.app.tools.checkv.checkv import CheckV
from camel.app.tools.checkv.checkvreporter import CheckVReporter
from camel.resources import DIR_CITATIONS


@dataclasses.dataclass(frozen=True)
class CheckVOpts(BaseOptions):
    """
    Custom options for CheckV.
    """
    working_dir: Path = dataclasses.field(default=Path.cwd(), metadata={'help': 'Working directory'})
    threads: int = dataclasses.field(default=4, metadata={'help': 'Number of threads to use'})

class MainCheckV(BaseScript[FastaInput, ScriptOutput, CheckVOpts]):
    """
    This class contains the main script for the CheckV tool.
    """

    def __init__(self, script_in: FastaInput, script_out: ScriptOutput, opts: CheckVOpts) -> None:
        """
        Initializes the main script.
        :param script_in: Script input
        :param script_out: Script output
        :param opts: Script options
        :return: None
        """
        super().__init__(
            name='CheckV',
            version='1.0',
            script_in=script_in,
            script_out=script_out,
            script_opts=opts
        )

    def _execute(self) -> None:
        """
        Runs the tool.
        :return: None
        """
        # Initialize the report
        report = reportutils.init_report(
            path_out=self._script_out.html, key='CheckV', title='CheckV', dir_out=self._script_out.dir)
        report.add_html_object(reportutils.create_overview_section(
            version=self.version, dataset_name=self._script_in.name, input_file_str=self._script_in.input_str))
        report.save()

        # Run CheckV
        script_in = self._script_in.create_symlinks(self._script_opts.working_dir / 'input')
        checkv = CheckV()
        checkv.add_input_files({'FASTA': [ToolIOFile(script_in.fasta)]})
        checkv.update_parameters(threads=self._script_opts.threads)
        checkv.run(self._script_opts.working_dir)

        # Create the output report
        checkv_reporter = CheckVReporter()
        checkv_reporter.add_input_files({key: checkv.tool_outputs[key] for key in checkv.tool_outputs.keys()})
        checkv_reporter.run(self._script_opts.working_dir)
        section = checkv_reporter.tool_outputs['HTML'][0].value
        report.add_html_object(section)
        section.copy_files(report.output_dir)

        # Add citation and command
        report.add_html_object(reportutils.create_commands_section(
            [checkv.informs], self._script_opts.working_dir))
        report.add_html_object(reportutils.create_citations_section(
            dir_=DIR_CITATIONS, keys_other=['Nayfach_2021-checkv']))
        report.save()


@click.command(name='checkv', short_help='Wrapper for CheckV')
@cliutils.add_click_options_from_dataclass(FastaInput)
@basescriptutils.add_output_opts
@cliutils.add_click_options_from_dataclass(CheckVOpts)
def main(**kwargs) -> None:
    """
    Runs the main script.
    :param kwargs: Keywords arguments
    :return: None
    """
    script = MainCheckV(
        script_in=FastaInput(**cliutils.from_kwargs(FastaInput, kwargs)),
        script_out=basescriptutils.parse_script_output(kwargs),
        opts=CheckVOpts(**cliutils.from_kwargs(CheckVOpts, kwargs))
    )
    script.run()

if __name__ == '__main__':
    initialize_logging()
    main()
