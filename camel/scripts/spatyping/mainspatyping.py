#!/usr/bin/env python
import argparse
import json
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.reports.htmlreport import HtmlReport
from camel.app.scriptutils import mainscriptutils
from camel.app.loggers import logger, initialize_logging
from camel.app.tools.blast.blastn import Blastn
from camel.app.tools.spatyping.spatyping import SpaTyping
from camel.app.tools.spatyping.spatypingreporter import SpaTypingReporter
from camel.app.wrappers.inputtype import helper_by_input_type


class MainSpaTyping:
    """
    This tool is used to run the Spa typing tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        :param args: Arguments, if not set they are removed from the command line
        """
        self._args = MainSpaTyping._parse_arguments(args)
        self._sample_name = mainscriptutils.determine_sample_name(self._args)
        self._helper = helper_by_input_type[self._args.input_type](self._args.working_dir, self._sample_name)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        mainscriptutils.add_common_arguments(argument_parser)
        mainscriptutils.add_input_files_arguments(argument_parser)
        mainscriptutils.add_assembly_arguments(argument_parser)
        argument_parser.add_argument(
            '--db-path', type=Path, help="Path to the database",
            default=Path('/db/pipelines/saureus/spatyping'))
        argument_parser.add_argument('--output-json', type=Path, help='Export output in JSON format')
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Initialize report
        report = mainscriptutils.init_report(
            self._args.output_html, self._args.output_dir, 'spa typing report', '<i>spa</i> typing')
        report.add_html_object(mainscriptutils.generate_analysis_info_section(self._args))
        report.save()

        # Run the tools
        fasta_file = self._helper.prepare_fasta_input(report, self._args)
        blastn_tsv_output = self.__run_blastn(fasta_file)
        spa_typing = self.__run_spa_tying(blastn_tsv_output.path, fasta_file)

        # Save JSON output (if specified)
        if self._args.output_json is not None:
            with open(self._args.output_json, 'w') as handle:
                json.dump(spa_typing.informs, handle, indent=2)
            logger.info(f'spa typing informs output save to: {self._args.output_json}')

        self.__add_report_output(spa_typing, report)

    def __run_blastn(self, fasta_file: Path) -> ToolIOFile:
        """
        Runs the BLASTN alignment.
        :param fasta_file: Input FASTA file
        :return: None
        """
        blastn = Blastn()
        blastn.add_input_files({
            'DB_BLAST': [ToolIOFile(self._args.db_path / 'profiles.fasta')],
            'FASTA': [ToolIOFile(fasta_file)]})
        blastn.update_parameters(
            output_format=SpaTyping.BLASTN_OUTPUT_FORMAT,
            num_alignments=100_000,
            task='blastn',
            dust='no'
        )
        blastn.run(self._args.working_dir)
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
            'CSV_profiles': [ToolIOFile(self._args.db_path / 'spatypes.csv')]
        })
        spa_typing.run(self._args.working_dir)
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
        reporter.run(self._args.working_dir)
        self._helper.export_output_and_commands_section(report, reporter.tool_outputs['VAL_HTML'][0].value)
        report.save()


def main(args: Sequence[str] | None = None) -> None:
    """
    Entry point for the common interface.
    :param args: Command line arguments
    :return: None
    """
    script = MainSpaTyping(args)
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
