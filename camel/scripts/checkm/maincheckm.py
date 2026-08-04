#!/usr/bin/env python
import dataclasses
import json
from datetime import datetime
from pathlib import Path

import click
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import reportutils

from camel.app.cli import cliutils
from camel.app.loggers import initialize_logging, logger
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.basescript.fastainput import FastaInput
from camel.app.scriptutils.model import BaseOptions, BaseOutput
from camel.app.tools.checkm.checkm import CheckM
from camel.app.tools.checkm.checkmreporter import CheckMReporter
from camel.resources import DIR_CITATIONS
from camel.version import __VERSION__


@dataclasses.dataclass(frozen=True)
class Output(BaseOutput):
    """
    Output for the script.
    """
    output_html: Path = dataclasses.field(metadata={'help': 'Output HTML file'})
    output_dir: Path | None = dataclasses.field(metadata={'help': 'Output directory'})
    output_json: Path | None = dataclasses.field(default=None, metadata={'help': 'Output TSV file'})


@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Custom options for the script.
    """
    working_dir: Path = dataclasses.field(default=Path.cwd(), metadata={'help': 'Working directory'})
    reduced_tree: bool = dataclasses.field(default=False, metadata={'help': 'Use reduced tree'})
    threads: int | None = dataclasses.field(default=4, metadata={'help': 'Number of threads to use'})


class MainCheckM(BaseScript[FastaInput, Output, Options]):
    """
    Main script for  the CheckM tool.
    """

    def __init__(self, in_: FastaInput, out: Output, opts: Options) -> None:
        """
        Initializes the main script.
        :param in_: Script input
        :param out: Script output
        :param opts: Options
        :return: None
        """
        tool_version = CheckM().version
        super().__init__(
            name='CheckM+',
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
            path_out=self._script_out.output_html,
            dir_out=self._script_out.output_dir,
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

        # Run CheckM
        checkm = CheckM()
        checkm.add_input_files(self.__prepare_input())
        checkm.update_parameters(
            reduced_tree=self._script_opts.reduced_tree, threads=self._script_opts.threads)
        checkm.run(self._script_opts.working_dir / 'checkm')

        # Save informs (if specified)
        if self._script_out.output_json is not None:
            with self._script_out.output_json.open('w') as handle:
                json.dump(checkm.informs, handle, indent=2)
                logger.info(f'CheckM informs saved to {self._script_out.output_json}')

        # Create the output report
        checkm_reporter = CheckMReporter()
        checkm_reporter.add_input_informs({'checkm': checkm.informs})
        checkm_reporter.add_input_files({'TSV': checkm.tool_outputs['TSV']})
        checkm_reporter.run(self._script_opts.working_dir)
        section = checkm_reporter.tool_outputs['HTML'][0].value
        section.copy_files(report.output_dir)
        report.add_html_object(section)

        # Add citation and command
        report.add_html_object(reportutils.create_commands_section(
            [checkm.informs], self._script_opts.working_dir))
        report.add_html_object(reportutils.create_citations_section(
            dir_=DIR_CITATIONS,
            keys_other=['Parks_2015-checkm']
        ))
        report.save()

    def __prepare_input(self) -> dict[str, list[ToolIOFile]]:
        """
        Prepares the input for the CheckM tool.
        :return: Input dictionary
        """
        dir_symlink = self._script_opts.working_dir / 'input'
        dir_symlink.mkdir(parents=True, exist_ok=True)
        script_in = self._script_in.create_symlinks(dir_symlink)
        input_dict = {'FASTA': [ToolIOFile(script_in.fasta)]}
        return input_dict


@click.command(name='checkm', short_help='Wrapper for CheckM')
@cliutils.add_click_options_from_dataclass(FastaInput)
@cliutils.add_click_options_from_dataclass(Output)
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Wrapper for CheckM.
    """
    script = MainCheckM(
        in_=FastaInput(**cliutils.from_kwargs(FastaInput, kwargs)),
        out=Output(**cliutils.from_kwargs(Output, kwargs)),
        opts=Options(**cliutils.from_kwargs(Options, kwargs))
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
