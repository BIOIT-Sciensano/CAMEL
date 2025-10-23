#!/usr/bin/env python
import argparse
import json
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional

from camel.app.core.reports import reportutils
from camel.app.core.reports.htmlreport import HtmlReport
from camel.app.scriptutils.basescript import BaseScript
from camel.app.scriptutils import mainscriptutils, absolute_path_by_pathlib
from camel.app.loggers import initialize_logging
from camel.app.wrappers.genedetectionwrapper import (
    GeneDetectionWrapper,
    GeneDetectionOutput,
)
from camel.app.wrappers.inputtype import helper_by_input_type


class MainGeneDetection(BaseScript):
    """
    This class is used to run the gene detection tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        :param args: Arguments (optional)
        """
        super().__init__(name='gene detection', version='1.0', snakefile=None)
        self._args = MainGeneDetection._parse_arguments(args)
        self._sample_name = mainscriptutils.determine_sample_name(self._args)
        self._helper = helper_by_input_type[self._args.input_type](self._args.working_dir, self._sample_name)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        mainscriptutils.add_common_arguments(argument_parser)
        mainscriptutils.add_assembly_arguments(argument_parser)
        mainscriptutils.add_input_files_arguments(argument_parser)
        argument_parser.add_argument('--database-dir', type=absolute_path_by_pathlib, required=True)
        argument_parser.add_argument('--detection-method', type=str, choices=['blast', 'kma'], default='blast')

        # BLAST specific parameters
        argument_parser.add_argument('--output-fasta', type=absolute_path_by_pathlib, help='output path for assembled contigs')
        argument_parser.add_argument('--blast-min-percent-identity', type=int, default=90)
        argument_parser.add_argument('--blast-min-percent-coverage', type=int, default=60)
        argument_parser.add_argument('--blast-task', type=str, choices=['blastn', 'megablast'], default='megablast')
        argument_parser.add_argument('--blast-filtering-method', type=str, choices=['cluster', 'score', 'overlap'], default='cluster')
        argument_parser.add_argument('--blast-score-nb-of-hits', type=int, default=5)
        argument_parser.add_argument('--blast-reads', action='store_true', default=None,
                                     help='perform blast search of the reads directly instead of on the assembly' )

        # KMA specific parameters
        argument_parser.add_argument('--kma-min-percent-identity', type=int, default=90)
        argument_parser.add_argument('--kma-min-percent-coverage', type=int, default=60)
        argument_parser.add_argument('--kma-ont', action='store_true', default=None)
        argument_parser.add_argument('--kma-cge', action='store_true', default=None)
        argument_parser.add_argument('--kma-apm', type=str, choices=['u', 'p', 'f'], default=None)
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        mainscriptutils.validate_input_files(self._args)

        # Initialize report
        report = reportutils.init_report(
            path_out=Path(self._args.output_html),
            key='Gene detection',
            title='Gene detection',
            dir_out=self._args.output_dir
        )
        report.add_html_object(reportutils.create_overview_section(
            version=self._version,
            dataset_name=self._sample_name,
            input_file_str=mainscriptutils.determine_input_file_str(self._args),
            extra_data=[('Detection method', self._args.detection_method)]
        ))
        report.save()

        # Prepare wrapper
        wrapper = GeneDetectionWrapper(self._helper.working_dir)
        db_data = self.__get_db_metadata()

        # Run wrapper
        if self._args.detection_method == 'blast':
            if self._args.blast_reads:
                fasta_input = self._helper.prepare_fasta_read_input(report, self._args)
            else:
                fasta_input = self._helper.prepare_fasta_input(report, self._args)
            # Save assembly if specified
            if self._args.output_fasta is not None:
                shutil.copyfile(str(fasta_input), self._args.output_fasta)
            wrapper.run_blast(fasta_input, self._sample_name, db_data, self._args.threads)
        elif self._args.detection_method == 'kma':
            fastq_input = self._helper.prepare_fastq_input(report, self._args)
            wrapper.run_kma(fastq_input, self._sample_name, db_data, self._args.threads)

        # Export all output
        self.__export_output(report, wrapper.output)

    def __get_db_metadata(self) -> dict[str, Any]:
        """
        Returns the database information dictionary.
        :return: Database information dictionary
        """
        config_data = {'path': str(self._args.database_dir)}

        # Add specific options
        if self._args.detection_method == 'blast':
            mainscriptutils.dict_merge(
                config_data,
                {
                    'params': {
                        'blastn': {
                            'blast_reads': True if self._args.blast_reads else False,
                            'filtering_method': self._args.blast_filtering_method,
                            'min_coverage': self._args.blast_min_percent_coverage,
                            'min_percent_identity': self._args.blast_min_percent_identity,
                            'score_nb_of_hits': self._args.blast_score_nb_of_hits,
                            'task': self._args.blast_task
                        }
                    }
                }
            )
        elif self._args.detection_method == 'kma':
            mainscriptutils.dict_merge(
                config_data,
                {
                    'params': {
                        'kma': {
                            'min_percent_identity': self._args.kma_min_percent_identity,
                            'min_coverage': self._args.kma_min_percent_coverage,
                            'ont': self._args.kma_ont,
                            'cge': self._args.kma_cge,
                            'apm': self._args.kma_apm
                        }
                    }
                }
            )

        # Add the extra metadata column
        with (self._args.database_dir / 'db_metadata.txt').open() as handle:
            db_metadata = json.load(handle)
            if 'extra_column' in db_metadata:
                config_data['metadata'] = db_metadata['extra_column']
        return config_data

    def __export_output(self, report: HtmlReport, output: GeneDetectionOutput) -> None:
        """
        Exports the output of the workflow.
        :param report: HTML report
        :param output: Workflow output
        :return: None
        """
        self._helper.logs['gene_detection'] = str(output.log_file) if output.log_file is not None else None
        self._helper.informs.append(output.informs)
        self._helper.export_output_and_commands_section(report, output.report_section)


if __name__ == '__main__':
    initialize_logging()
    main = MainGeneDetection()
    main.run()
