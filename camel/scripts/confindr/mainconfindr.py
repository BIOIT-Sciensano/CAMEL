#!/usr/bin/env python
import argparse
import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional

from camel.app.core.reports import reportutils
from camel.app.scriptutils.basescript import BaseScript
from camel.app.scriptutils import mainscriptutils
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger, initialize_logging
from camel.app.tools.confindr.confindr import ConFindr
from camel.app.tools.confindr.confindrreporter import ConFindrReporter


class MainConFindr(BaseScript):
    """
    This class contains the main script for the ConFindr tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        """
        super().__init__(name='ConFindr', version='1.0', snakefile=None)
        self._args = MainConFindr._parse_arguments(args)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        mainscriptutils.add_input_files_arguments(argument_parser, False)
        mainscriptutils.add_common_arguments(argument_parser)
        argument_parser.add_argument('--db', type=Path, required=True)
        argument_parser.add_argument(
            '--output-json', type=Path, help='If specified, ConFindr informs are stored in this file')

        # Parameters
        argument_parser.add_argument(
            '--rmlst', action='store_true', help='Prefer using rMLST databases over core-gene derived databases')
        argument_parser.add_argument('--quality-cutoff', type=int, default=20, help='Base quality cutoff')
        argument_parser.add_argument('--base-cutoff', type=int, default=3, help='Number of bases  cutoff')
        argument_parser.add_argument('--base-percentage-cutoff', type=int, default=5, help='Base percentage cutoff')
        argument_parser.add_argument(
            '--min-matching-hashes', type=int, default=150, help='Minimum number of matching KMA hashes')

        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the tool.
        :return: None
        """
        # Initialize the report
        report = reportutils.init_report(
            path_out=self._args.output_html,
            key=self._name,
            title='ConFindr',
            dir_out=self._args.output_dir)
        report.add_html_object(reportutils.create_overview_section(
            version=self.version,
            dataset_name=self._name,
            input_file_str=mainscriptutils.determine_input_file_str(self._args),
            input_type=self._args.input_type
        ))
        report.save()

        # Check if the database exists
        if not self._args.db.exists():
            raise FileNotFoundError(f'DB not found: {self._args.db}')

        # Run ConFindr
        input_dict = self.__prepare_input()
        confindr = ConFindr()
        confindr.update_parameters(
            databases=str(self._args.db),
            quality_cutoff=self._args.quality_cutoff,
            base_cutoff=self._args.base_cutoff,
            base_fraction_cutoff=self._args.base_percentage_cutoff / 100,
            min_matching_hashes=self._args.min_matching_hashes,
            data_type={'ont': 'Nanopore', 'illumina': 'Illumina'}[self._args.input_type],
            threads=self._args.threads
        )
        if self._args.rmlst:
            confindr.update_parameters(rmlst=True)
        confindr.add_input_files(input_dict)
        confindr.run(self._args.working_dir)

        # Save informs to file (if specified)
        if self._args.output_json is not None:
            with self._args.output_json.open('w') as handle:
                json.dump(confindr.informs, handle, indent=2)
                logger.info(f'ConFindr informs saved to {self._args.output_json}')

        # Create output report
        confindr_reporter = ConFindrReporter()
        confindr_reporter.add_input_informs({'confindr': confindr.informs})
        confindr_reporter.update_parameters(input_type=self._args.input_type)
        confindr_reporter.run(self._args.working_dir)
        report.add_html_object(confindr_reporter.tool_outputs['HTML'][0].value)

        # Add citation and command
        report.add_html_object(reportutils.create_commands_section([confindr.informs], self._args.working_dir))
        report.add_html_object(reportutils.create_citations_section(['Low_2019-confindr', 'Jolley_2012-rmlst']))
        report.save()
        logger.info(f'Report saved to: {self._args.output_html}')

    def __prepare_input(self) -> dict[str, Any]:
        """
        Prepares the input for the confindr tool.
        :return: Input dictionary
        """
        if self._args.fastq_pe is not None:
            return {'FASTQ_PE': [ToolIOFile(fq) for fq in self._args.fastq_pe]}
        else:
            return {'FASTQ_SE': [ToolIOFile(self._args.fastq_se)]}


if __name__ == '__main__':
    initialize_logging()
    main = MainConFindr()
    main.run()
