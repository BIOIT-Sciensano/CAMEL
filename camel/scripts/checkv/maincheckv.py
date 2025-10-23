#!/usr/bin/env python
import argparse
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional

from camel.app.core.reports import reportutils
from camel.app.core.utils import fileutils
from camel.app.scriptutils import mainscriptutils
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import initialize_logging
from camel.app.tools.checkv.checkv import CheckV
from camel.app.tools.checkv.checkvreporter import CheckVReporter


class MainCheckV:
    """
    This class contains the main script for the CheckV tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        """
        self._args = MainCheckV._parse_arguments(args)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--fasta', help='FASTA input', type=Path, required=True)
        argument_parser.add_argument('--fasta-name', help='FASTA input name')
        argument_parser.add_argument('--working-dir', help='Working directory', type=Path, default=Path.cwd())
        argument_parser.add_argument('--output-html', help='Report output', type=Path)
        argument_parser.add_argument('--output-dir', help='Output directory', type=Path)
        argument_parser.add_argument('--threads', type=int, default=4, help='Number of threads to use')
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the tool.
        :return: None
        """
        input_dict, fasta_name = self.__prepare_input()
        report = mainscriptutils.init_report(self._args.output_html, self._args.output_dir, 'CheckV', 'CheckV')
        report.add_html_object(mainscriptutils.generate_analysis_info_section(self._args, input_file_str=fasta_name))
        report.save()

        # Run CheckV
        checkv = CheckV()
        checkv.add_input_files(input_dict)
        checkv.update_parameters(threads=self._args.threads)
        checkv.run(self._args.working_dir)

        # Create output report
        checkv_reporter = CheckVReporter()
        checkv_reporter.add_input_files({key: checkv.tool_outputs[key] for key in checkv.tool_outputs.keys()})
        checkv_reporter.run(self._args.working_dir)
        section = checkv_reporter.tool_outputs['HTML'][0].value
        report.add_html_object(section)
        section.copy_files(report.output_dir)

        # Add citation and command
        report.add_html_object(reportutils.create_commands_section([checkv.informs], self._args.working_dir))
        report.add_html_object(reportutils.create_citations_section(['Parks_2015-checkm']))
        report.save()

    def __prepare_input(self) -> tuple[dict[str, Any], str]:
        """
        Prepares the input for the CheckM tool.
        :return: Input dictionary
        """
        dir_input = self._args.working_dir / 'input'
        dir_input.mkdir(parents=True, exist_ok=True)
        fasta_name = self._args.fasta_name if self._args.fasta_name else self._args.fasta.name
        path_new = dir_input / fileutils.make_valid(fasta_name)
        path_new.symlink_to(self._args.fasta)
        return {'FASTA': [ToolIOFile(path_new)]}, fasta_name


if __name__ == '__main__':
    initialize_logging()
    main = MainCheckV()
    main.run()
