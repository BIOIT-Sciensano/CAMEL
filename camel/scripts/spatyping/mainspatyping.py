#!/usr/bin/env python
import dataclasses
import json
from pathlib import Path

import click
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.reports.htmlreport import HtmlReport
from camelcore.app.utils import reportutils

from camel.app.cli import cliutils
from camel.app.loggers import initialize_logging, logger
from camel.app.scriptutils import inputhelper, model
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.app.scriptutils.inputhelper import helper_by_input_type
from camel.app.scriptutils.inputhelper.inputhelperbase import InputHelperBase
from camel.app.scriptutils.model import BaseOptions
from camel.app.tools.blast.blastn import Blastn
from camel.app.tools.spatyping.spatyping import SpaTyping
from camel.app.tools.spatyping.spatypingreporter import SpaTypingReporter
from camel.version import __VERSION__


@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Options for the spa typing script.
    """
    db: Path = dataclasses.field(metadata={'help': 'Path to the database'})
    working_dir: Path = dataclasses.field(default=Path('working'), metadata={
        'help': 'Working directory', 'show_default': True})
    threads: int | None = dataclasses.field(default=1, metadata={'help': 'Number of threads', 'show_default': True})


class MainSpaTyping(BaseScript[ScriptInput, ScriptOutput, Options]):
    """
    Main script for the spa typing tool.
    """

    def __init__(
            self, script_in: ScriptInput, helper: InputHelperBase, script_out: ScriptOutput,
            script_opts: Options) -> None:
        """
        Initializes the main script.
        :param script_in: Script input
        :param helper: Input helper
        :param script_out: Script output
        :param script_opts: Script options
        """
        super().__init__(
            name='spa typing',
            title='<i>spa</i> typing',
            version=f'CAMEL_{__VERSION__}',
            script_in=script_in,
            script_out=script_out,
            script_opts=script_opts
        )
        self._helper: InputHelperBase = helper

    def _execute(self) -> None:
        """
        Executes the script.
        :return: None
        """
        # Initialize report
        report = reportutils.init_report(
            path_out=Path(self._script_out.html),
            key=self.name,
            title=self.title,
            dir_out=self._script_out.dir
        )
        report.add_html_object(reportutils.create_overview_section(
            version=self._version,
            dataset_name=self._script_in.name,
            input_file_str=self._script_in.input_str
        ))
        report.save()

        # Run the tools
        fasta_file = self._helper.prepare_fasta_input(self._script_in, report)
        blastn_tsv_output = self.__run_blastn(fasta_file)
        spa_typing = self.__run_spa_tying(blastn_tsv_output.path, fasta_file)

        # Save JSON output (if specified)
        if self._script_out.json is not None:
            with open(self._script_out.json, 'w') as handle:
                json.dump(spa_typing.informs, handle, indent=2)
            logger.info(f'spa typing informs output save to: {self._script_out.json}')
        self.__add_report_output(spa_typing, report)

    def __run_blastn(self, fasta_file: Path) -> ToolIOFile:
        """
        Runs the BLASTN alignment.
        :param fasta_file: Input FASTA file
        :return: None
        """
        blastn = Blastn()
        blastn.add_input_files({
            'DB_BLAST': [ToolIOFile(self._script_opts.db / 'profiles.fasta')],
            'FASTA': [ToolIOFile(fasta_file)]})
        blastn.update_parameters(
            dust='no',
            num_alignments=100_000,
            output_format=SpaTyping.BLASTN_OUTPUT_FORMAT,
            task='blastn',
            threads=self._script_opts.threads,
        )
        blastn.run(self._script_opts.working_dir)
        self._helper.informs.append(blastn.informs)
        return blastn.tool_outputs['TSV'][0]

    def __run_spa_tying(self, tsv_output: Path, fasta_file: Path) -> SpaTyping:
        """
        Runs the Spa typing on the tabular blast output.
        :param tsv_output: Tabular blast output
        :param fasta_file: FASTA file
        :return: None
        """
        spa_typing = SpaTyping()
        spa_typing.add_input_files({
            'TSV': [ToolIOFile(tsv_output)],
            'FASTA': [ToolIOFile(fasta_file)],
            'CSV_profiles': [ToolIOFile(self._script_opts.db / 'spatypes.csv')]
        })
        spa_typing.run(self._script_opts.working_dir)
        return spa_typing

    def __add_report_output(self, spa_typing: SpaTyping, report: HtmlReport) -> None:
        """
        Adds the spa typing output to the report.
        :param spa_typing: Spa typing tool instance
        :param report: Report to append information to
        :return: None
        """
        reporter = SpaTypingReporter()
        reporter.add_input_informs({'spa_typing': spa_typing.informs})
        reporter.add_input_files({'VAL_hits': spa_typing.tool_outputs['VAL_hits']})
        reporter.run(self._script_opts.working_dir)
        self._helper.export_output_and_commands_section(report, reporter.tool_outputs['VAL_HTML'][0].value)
        report.save()

@click.command(name='spa_typing', short_help='spa typing for Staphylococcus aureus')
@basescriptutils.add_input_opts(supported=[model.InputType.FASTA, model.InputType.ILLUMINA, model.InputType.ONT])
@basescriptutils.add_output_opts
@cliutils.add_click_options_from_dataclass(Options)
@inputhelper.add_helper_opts
def main(**kwargs) -> None:
    """
    Spa typing for Staphylococcus aureus.
    """
    # Parse the script input
    script_input = basescriptutils.parse_script_input(kwargs)
    script_opts = Options(**cliutils.from_kwargs(Options, kwargs))

    # Initialize the helper class to prepare the input
    helper = helper_by_input_type[script_input.type_](dir_=script_opts.working_dir, name=script_input.name)
    helper.set_opts(*helper.opts_from_cli(kwargs))

    # Run the main script
    script = MainSpaTyping(
        script_in=basescriptutils.parse_script_input(kwargs),
        script_out=basescriptutils.parse_script_output(kwargs),
        script_opts=script_opts,
        helper=helper
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
