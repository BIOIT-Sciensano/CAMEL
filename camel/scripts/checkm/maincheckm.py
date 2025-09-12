#!/usr/bin/env python
import argparse
import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.checkm.checkm import CheckM
from camel.app.tools.checkm.checkmreporter import CheckMReporter


class MainCheckM:
    """
    This class contains the main script for the CheckM tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        """
        self._args = MainCheckM._parse_arguments(args)
        self._camel = Camel()

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--fasta', nargs=2, action='append', help='FASTA input', required=True)
        argument_parser.add_argument('--working-dir', help='Working directory', type=Path, default=Path.cwd())
        argument_parser.add_argument('--output-html', type=Path, help='Report output')
        argument_parser.add_argument('--output-dir', type=Path, help='Output directory')
        argument_parser.add_argument('--output-json', type=Path, help='Output path to store CheckM informs')
        argument_parser.add_argument('--threads', type=int, default=4, help='Number of threads to use')
        argument_parser.add_argument('--reduced_tree', action='store_true', default=None)
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the tool.
        :return: None
        """
        input_dict, input_files_str = self.__prepare_input()

        # Initialize report
        report = mainscriptutils.init_report(self._args.output_html, self._args.output_dir, 'CheckM', 'CheckM')
        report.add_html_object(mainscriptutils.generate_analysis_info_section(
            self._args, input_file_str=input_files_str))
        report.save()

        # Run CheckM
        checkm = CheckM()
        checkm.add_input_files(input_dict)
        checkm.update_parameters(threads=self._args.threads)
        checkm.update_parameters(reduced_tree=self._args.reduced_tree)
        checkm.run(self._args.working_dir)

        # Save informs (if specified)
        if self._args.output_json is not None:
            with self._args.output_json.open('w') as handle:
                json.dump(checkm.informs, handle, indent=2)
                logger.info(f'CheckM informs saved to {self._args.output_json}')

        # Create output report
        checkm_reporter = CheckMReporter()
        checkm_reporter.add_input_informs({'checkm': checkm.informs})
        checkm_reporter.add_input_files({'TSV': checkm.tool_outputs['TSV']})
        checkm_reporter.run(self._args.working_dir)
        section = checkm_reporter.tool_outputs['HTML'][0].value
        section.copy_files(report.output_dir)
        report.add_html_object(section)

        # Add citation and command
        report.add_html_object(SnakePipelineUtils.create_commands_section([checkm.informs], self._args.working_dir))
        report.add_html_object(SnakePipelineUtils.create_citations_section(['Parks_2015-checkm']))
        report.save()

    def __prepare_input(self) -> tuple[dict[str, Any], str]:
        """
        Prepares the input for the CheckM tool.
        :return: Input dictionary
        """
        self._args.working_dir.mkdir(parents=True, exist_ok=True)
        input_dict = {'FASTA': []}
        for fasta_file, fasta_name in self._args.fasta:
            path_new = self._args.working_dir / FileSystemHelper.make_valid(fasta_name)
            path_new.symlink_to(fasta_file)
            input_dict['FASTA'].append(ToolIOFile(path_new))
        input_str = ', '.join([fasta_name for _, fasta_name in self._args.fasta])
        return input_dict, input_str


if __name__ == '__main__':
    Camel.get_instance()
    main = MainCheckM()
    main.run()
