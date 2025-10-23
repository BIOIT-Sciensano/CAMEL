#!/usr/bin/env python
import argparse
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

from camel.app.core.reports import reportutils
from camel.app.core.utils import fileutils
from camel.app.scriptutils import mainscriptutils
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import initialize_logging
from camel.app.tools.integronfinder.integronfinder import IntegronFinder
from camel.app.tools.integronfinder.integronfinderreporter import IntegronFinderReporter


class MainIntegronFinder:
    """
    Main script for the IntegronFinder tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        :param args: Arguments (optional)
        """
        self._args = MainIntegronFinder._parse_arguments(args)
        self._sample_name = mainscriptutils.determine_sample_name(self._args)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments to parse
        :return: Parsed arguments
        """
        parser = argparse.ArgumentParser()
        mainscriptutils.add_common_arguments(parser)
        parser.add_argument('--fasta', help='Input FASTA file', type=Path, required=True)
        parser.add_argument('--fasta-name', help="Input FASTA file name (For Galaxy)", type=str)
        parser.add_argument('--output-tsv', help="Copy the tabular file to this location", type=Path)
        parser.add_argument('--local-max', action='store_true',
                            help='Allows thorough local detection (slower but more sensitive).')
        return parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        # Initialize report
        report = mainscriptutils.init_report(
            self._args.output_html, self._args.output_dir, 'IntegronFinder report', 'IntegronFinder')
        input_file_str = self._args.fasta_name if self._args.fasta_name is not None else self._args.fasta.name
        report.add_html_object(mainscriptutils.generate_analysis_info_section(
            self._args, input_file_str=input_file_str))
        report.save()

        # Run IntegronFinder
        integron_finder = IntegronFinder()
        integron_finder.add_input_files({'FASTA': [ToolIOFile(self._args.fasta)]})
        integron_finder.update_parameters(threads=self._args.threads)
        if self._args.local_max is True:
            integron_finder.update_parameters(local_max=True)
        integron_finder.run(self._args.working_dir)

        # Create output section
        reporter = IntegronFinderReporter()
        reporter.add_input_files({'TSV': integron_finder.tool_outputs['TSV']})
        reporter.add_input_informs({'integron_finder': integron_finder.informs})
        reporter.update_parameters(name=fileutils.make_valid(self._sample_name))
        reporter.run(self._args.working_dir)
        report.add_html_object(reporter.tool_outputs['HTML'][0].value)
        reporter.tool_outputs['HTML'][0].value.copy_files(report.output_dir)

        # Add commands and citations
        report.add_html_object(reportutils.create_commands_section(
            [integron_finder.informs], self._args.working_dir))
        report.add_html_object(reportutils.create_citations_section(['Neron_2022-integronfinder']))
        report.save()

        # Copy the TSV output file when specified
        if self._args.output_tsv is not None:
            shutil.copyfile(integron_finder.tool_outputs['TSV'][0].path, self._args.output_tsv)


if __name__ == '__main__':
    initialize_logging()
    main = MainIntegronFinder()
    main.run()
