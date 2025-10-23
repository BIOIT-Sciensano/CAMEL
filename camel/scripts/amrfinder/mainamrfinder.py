#!/usr/bin/env python
import argparse
import shutil
from collections.abc import Sequence
from datetime import datetime
from pathlib import Path
from typing import Optional

from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.reports import reportutils
from camel.app.scriptutils.basescript import BaseScript
from camel.app.core.utils import fileutils
from camel.app.scriptutils import mainscriptutils
from camel.app.loggers import initialize_logging
from camel.app.tools.amrfinder.amrfinder import AMRFinder
from camel.app.tools.amrfinder.amrfinderreporter import AMRFinderReporter


class MainAMRFinder(BaseScript):
    """
    Main script for the AMRFinder tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        :param args: Arguments (optional)
        """
        super().__init__(name='AMRFinder+', version='1.0', snakefile=None)
        self._args = MainAMRFinder._parse_arguments(args)
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
        parser.add_argument('--db', help='Database', type=Path, required=True)
        parser.add_argument('--min-cov', help='Minimum target coverage', type=int, default=50)
        parser.add_argument(
            '--min-id', type=int, help='Minimum % identity (defaults to curated cutoffs when not specified)')
        parser.add_argument(
            '--organism', help='Organisms, required to identify point mutations associated with AMR', type=str)
        return parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        # Initialize report
        report = reportutils.init_report(
            path_out=self._args.output_html,
            dir_out=self._args.output_dir,
            key=self.name,
            title=f'{self.name} v{self.version}',
        )

        # Create the overview section
        additional_info = [
            ['Organism', self._args.organism if self._args.organism is not None else 'Not specified'],
            ['Min % identity', self._args.min_id if self._args.min_id is not None else 'Curated (default)'],
            ['Min % coverage', self._args.min_cov],
        ]
        input_file_str = self._args.fasta_name if self._args.fasta_name is not None else self._args.fasta.name
        report.add_html_object(reportutils.create_overview_section(
            version=self.version,
            dataset_name=self._sample_name,
            input_file_str=input_file_str,
            date=datetime.now(),
            extra_data=additional_info
        ))
        report.save()

        # Run AMRFinder
        amrfinder = AMRFinder()
        amrfinder.add_input_files({
            'FASTA': [ToolIOFile(self._args.fasta)],
            'DIR': [ToolIODirectory(self._args.db)]
        })
        amrfinder.update_parameters(
            min_cov=self._args.min_cov / 100.0,
            output_path=f'amrfinder_{fileutils.make_valid(self._sample_name)}.tsv'
        )
        if self._args.min_id is not None:
            amrfinder.update_parameters(min_ident=self._args.min_id / 100.0)
        if self._args.organism is not None:
            amrfinder.update_parameters(organism=self._args.organism)
        self._args.working_dir.mkdir(parents=True, exist_ok=True)
        amrfinder.run(self._args.working_dir)

        # Create the output section
        amrfinder_reporter = AMRFinderReporter()
        amrfinder_reporter.add_input_files({'TSV': amrfinder.tool_outputs['TSV']})
        amrfinder_reporter.add_input_informs({'amrfinder': amrfinder.informs})
        amrfinder_reporter.run(self._args.working_dir)
        report.add_html_object(amrfinder_reporter.tool_outputs['HTML'][0].value)
        amrfinder_reporter.tool_outputs['HTML'][0].value.copy_files(report.output_dir)

        # Add commands and citations
        report.add_html_object(reportutils.create_commands_section([amrfinder.informs], self._args.working_dir))
        report.add_html_object(reportutils.create_citations_section(['Feldgarden_2019-ndaro']))
        report.save()

        # Copy the TSV output file when specified
        if self._args.output_tsv is not None:
            shutil.copyfile(amrfinder.tool_outputs['TSV'][0].path, self._args.output_tsv)


if __name__ == '__main__':
    initialize_logging()
    main = MainAMRFinder()
    main.run()
