#!/usr/bin/env python
import argparse
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

import pandas as pd

from camel.app.core.reports import reportutils
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.scriptutils import mainscriptutils
from camel.app.loggers import initialize_logging
from camel.app.tools.btyper.btyper import BTyper
from camel.app.tools.btyper.btyperreporter import BTyperReporter


class MainBTyper:
    """
    This class is used to run the main ResFinder local script.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        :param args: Arguments (optional)
        """
        self._args = MainBTyper.parse_arguments(args)
        self._sample_name = mainscriptutils.determine_sample_name(self._args)

    @staticmethod
    def parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        mainscriptutils.add_common_arguments(argument_parser)
        argument_parser.add_argument('--fasta', help='Input FASTA file', type=Path, required=True)
        argument_parser.add_argument('--fasta-name', help='Input FASTA file name', type=str)
        argument_parser.add_argument('--virulence', help='perform virulence gene detection', action='store_true')
        argument_parser.add_argument('--bt', help='perform Bt toxin gene detection for cry, cyt, and vip genes',
                                     action='store_true')
        argument_parser.add_argument('--mlst', help='assign genome to a sequence type', action='store_true')
        argument_parser.add_argument('--panc', help='assign genome to a phylogenetic group using an adjusted, '
                                                    'eight-group panC group assignment scheme', action='store_true')
        argument_parser.add_argument('--output-tsv', help='Copy the output tabular file to this location', type=Path)
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        # Initialize report
        report = mainscriptutils.init_report(
            self._args.output_html, self._args.output_dir, 'BTyper3 report', 'BTyper3')
        additional_info = [
            ['Virulence:', '{}'.format('True' if self._args.virulence else 'False')],
            ['MLST:', '{}'.format('True' if self._args.mlst else 'False')],
            ['PanC:', '{}'.format('True' if self._args.panc else 'False')],
            ['BT:', '{}'.format('True' if self._args.bt else 'False')],
        ]
        report.add_html_object(mainscriptutils.generate_analysis_info_section(self._args, additional_info))
        report.save()

        # Run tools
        btyper = self.__run_btyper()
        section = self.__run_reporter(btyper)
        report.add_html_object(section)
        section.copy_files(report.output_dir)

        # Save report
        all_informs = [btyper.informs]
        report.add_html_object(reportutils.create_commands_section(all_informs, self._args.working_dir))
        report.add_html_object(reportutils.create_citations_section([
            'Carroll_2020a-btyper3', 'Carroll_2020b-btyper3']))
        report.save()

        # Copy the TSV output file when specified
        if self._args.output_tsv is not None:
            if self._args.fasta_name:
                self.__add_fasta_name_to_tsv(btyper.tool_outputs['TSV'][0].path)
            shutil.copyfile(btyper.tool_outputs['TSV'][0].path, self._args.output_tsv)

    def __run_btyper(self) -> BTyper:
        """
        Runs BTyper.
        :return: BTyper tool instance.
        """
        btyper = BTyper()
        btyper.add_input_files({'FASTA': [ToolIOFile(self._args.fasta)]})

        # Update parameters
        if not self._args.virulence:
            btyper.update_parameters(virulence='False')
        if not self._args.bt:
            btyper.update_parameters(bt='False')
        if not self._args.panc:
            btyper.update_parameters(panc='False')
        if not self._args.mlst:
            btyper.update_parameters(mlst='False')
        btyper.update_parameters(output_dir=self._args.output_dir)

        # Run the tool
        btyper.run(self._args.working_dir)
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
        data_in['#filename'] = f'{self._args.fasta_name}'
        data_in.to_csv(input_file, sep='\t', header=True)


if __name__ == '__main__':
    initialize_logging()
    btyper_main = MainBTyper()
    btyper_main.run()
