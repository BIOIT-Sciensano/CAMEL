#!/usr/bin/env python
import dataclasses
import shutil
from datetime import datetime
from pathlib import Path

import click
import pandas as pd

from camel.app.cli import cliutils
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.reports import reportutils
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.loggers import initialize_logging
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.basescript.fastainput import FastaInput
from camel.app.scriptutils.model import BaseOutput, BaseOptions
from camel.app.tools.btyper.btyper import BTyper
from camel.app.tools.btyper.btyperreporter import BTyperReporter


@dataclasses.dataclass(frozen=True)
class Output(BaseOutput):
    """
    Output for the BTyper script.
    """
    output_html: Path = dataclasses.field(metadata={'help': 'Output HTML file'})
    output_dir: Path | None = dataclasses.field(metadata={'help': 'Output directory'})
    output_tsv: Path | None = dataclasses.field(default=None, metadata={'help': 'Output TSV file'})


@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Options for the BTYper script.
    """
    working_dir: Path = dataclasses.field(default=Path('working'), metadata={
        'help': 'Working directory', 'show_default': True})
    virulence: bool = dataclasses.field(default=False, metadata={'help': 'Perform virulence gene detection'})
    bt: bool = dataclasses.field(default=False, metadata={'help': 'Perform Bt toxin gene detection'})
    mlst: bool = dataclasses.field(default=False, metadata={'help': 'Assign genome to a sequence type'})
    panc: bool = dataclasses.field(default=False, metadata={'help': 'Assign genome to a phylogenetic group'})


class MainBTyper(BaseScript[FastaInput, Output, Options]):
    """
    Main script for the BTyper tool.
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
            name='BTyper',
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
            path_out=self._script_out.output_html,
            dir_out=self._script_out.output_dir,
            key=self.name,
            title=f'{self.name} v{self.version}',
        )
        additional_info = [
            ['Virulence:', str(self._script_opts.virulence)],
            ['MLST:', str(self._script_opts.mlst)],
            ['PanC:', str(self._script_opts.panc)],
            ['BT:', str(self._script_opts.bt)],
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
        btyper = self.__run_btyper()
        section = self.__run_reporter(btyper)
        report.add_html_object(section)
        section.copy_files(report.output_dir)

        # Save report
        all_informs = [btyper.informs]
        report.add_html_object(reportutils.create_commands_section(
            all_informs, self._script_opts.working_dir))
        report.add_html_object(reportutils.create_citations_section([
            'Carroll_2020a-btyper3', 'Carroll_2020b-btyper3']))
        report.save()

        # Copy the TSV output file when specified
        if self._script_out.output_tsv is not None:
            self.__add_fasta_name_to_tsv(btyper.tool_outputs['TSV'][0].path)
            shutil.copyfile(btyper.tool_outputs['TSV'][0].path, self._script_out.output_tsv)

    def __run_btyper(self) -> BTyper:
        """
        Runs BTyper.
        :return: BTyper tool instance.
        """
        btyper = BTyper()
        btyper.add_input_files({'FASTA': [ToolIOFile(self._script_in.fasta)]})

        # Update parameters
        if not self._script_opts.virulence:
            btyper.update_parameters(virulence='False')
        if not self._script_opts.bt:
            btyper.update_parameters(bt='False')
        if not self._script_opts.panc:
            btyper.update_parameters(panc='False')
        if not self._script_opts.mlst:
            btyper.update_parameters(mlst='False')
        btyper.update_parameters(output_dir=self._script_out.output_dir)

        # Run the tool
        btyper.run(self._script_opts.working_dir)
        return btyper

    def __run_reporter(self, btyper: BTyper) -> HtmlReportSection:
        """
        Runs the BTyper reporter.
        :param btyper: BTyper tool instance.
        :return: None.
        """
        reporter = BTyperReporter()
        reporter.add_input_files({'TSV': btyper.tool_outputs['TSV']})
        reporter.add_input_informs({'btyper': btyper.informs})
        reporter.run()
        return reporter.tool_outputs['HTML'][0].value

    def __add_fasta_name_to_tsv(self, input_file: Path) -> None:
        """
        Modify the filename entry by the fasta file name. Important for Galaxy.
        :param input_file: TSV output of BTyper
        :return: None.
        """
        data_in = pd.read_table(input_file)
        data_in['#filename'] = self._script_in.input_str
        data_in.to_csv(input_file, sep='\t', header=True)


@click.command(name='btyper', short_help='Wrapper for BTyper')
@cliutils.add_click_options_from_dataclass(FastaInput)
@cliutils.add_click_options_from_dataclass(Output)
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Wrapper for BTyper.
    """
    script = MainBTyper(
        in_=FastaInput(**cliutils.from_kwargs(FastaInput, kwargs)),
        out=Output(**cliutils.from_kwargs(Output, kwargs)),
        opts=Options(**cliutils.from_kwargs(Options, kwargs))
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
