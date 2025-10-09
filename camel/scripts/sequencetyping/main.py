#!/usr/bin/env python
import argparse
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.pipelines import absolute_path_by_pathlib
from camel.app.components.sequencetyping.sequencetypingutils import SequenceTypingUtils
from camel.app.components.workflows.inputtype import helper_by_input_type
from camel.app.components.workflows.sequencetypingwrapper import (
    SequenceTypingInput,
    SequenceTypingOutput,
    SequenceTypingWrapper,
)
from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.io.tooliofile import ToolIOFile


class MainSequenceTyping:
    """
    Class to run sequence typing tool, it supports both BLAST+ and SRST2 as detection methods for alleles.
    """

    def __init__(self, args: Optional[Sequence[str]] = None):
        """
        Initializes the main script.
        :param args: (Optional) arguments
        """
        self._args = MainSequenceTyping._parse_arguments(args)
        self._sample_name = mainscriptutils.determine_sample_name(self._args)
        self._helper = helper_by_input_type[self._args.input_type](self._args.working_dir, self._sample_name)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: (Optional) arguments
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        mainscriptutils.add_common_arguments(argument_parser)
        mainscriptutils.add_assembly_arguments(argument_parser)
        mainscriptutils.add_input_files_arguments(argument_parser)
        argument_parser.add_argument('--scheme-dir', required=True, type=absolute_path_by_pathlib)
        argument_parser.add_argument(
            '--detection-method', type=str, choices=['blast', 'kma', 'mist'], default='blast')
        argument_parser.add_argument(
            '--output-fasta', type=absolute_path_by_pathlib,
            help='output path for assembled contigs (only used for BLAST-based detection)')
        argument_parser.add_argument(
            '--output-tsv', type=absolute_path_by_pathlib,
            help='Output path for the tabular output file (does not work with mixed schemes)')
        argument_parser.add_argument('--blastn-task', type=str, choices=['blastn', 'megablast'], default='megablast')
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the workflow.
        :return: None
        """
        mainscriptutils.validate_input_files(self._args)

        # Initialize report
        report = mainscriptutils.init_report(
            self._args.output_html, self._args.output_dir, 'Sequence typing report',f'Sequence typing {self._args.detection_method}')
        report.add_html_object(mainscriptutils.generate_analysis_info_section(self._args))
        report.save()

        # Run script with wrapper
        db_data = SequenceTypingUtils.parse_scheme_metadata(self._args.scheme_dir)
        if self._args.detection_method in ('blast', 'mist'):
            fasta_file = self._helper.prepare_fasta_input(report, self._args)
            # Save assembly if specified
            if self._args.output_fasta is not None:
                shutil.copyfile(str(fasta_file), self._args.output_fasta)
            if self._args.detection_method == 'blast':
                output = self.__run_sequence_typing_blast(fasta_file, db_data['name'], self._args.scheme_dir)
            else:
                output = self.__run_sequence_typing_mist(fasta_file, db_data['name'], self._args.scheme_dir)
        elif self._args.detection_method == 'kma':
            fastq_input = self._helper.prepare_fastq_input(report, self._args)
            output = self.__run_sequence_typing_kma(fastq_input, db_data['name'], self._args.scheme_dir)
        else:
            raise ValueError(f"Invalid detection method: {self._args.detection_method}")
        self.__export_output(output, report)

    def __run_sequence_typing_blast(self, fasta_file: Path, db_key: str, db_path: Path) -> SequenceTypingOutput:
        """
        Runs the sequence typing workflow using BLAST.
        :param fasta_file: Input FASTA file
        :param db_key: Database key
        :param db_path: Database directory path
        :return: None
        """
        wrapper = SequenceTypingWrapper(self._args.working_dir)
        workflow_input = SequenceTypingInput(
            sample_name=self._sample_name,
            fasta=ToolIOFile(fasta_file),
            input_type='fasta',
            db_path=db_path,
            db_key=db_key
        )
        wrapper.run_workflow_blast(workflow_input, self._args.blastn_task, self._args.threads)
        return wrapper.output

    def __run_sequence_typing_mist(self, fasta_file: Path, db_key: str, db_path: Path) -> SequenceTypingOutput:
        """
        Runs the sequence typing workflow using MiST.
        :param fasta_file: Input FASTA file
        :param db_key: Database key
        :param db_path: Database directory path
        :return: None
        """
        wrapper = SequenceTypingWrapper(self._args.working_dir)
        workflow_input = SequenceTypingInput(
            sample_name=self._sample_name,
            fasta=ToolIOFile(fasta_file),
            input_type='fasta',
            db_path=db_path,
            db_key=db_key
        )
        wrapper.run_workflow_mist(workflow_input, threads=self._args.threads)
        return wrapper.output

    def __run_sequence_typing_kma(self, fastq_input: FastqInput, db_key: str, db_path: Path) -> \
            SequenceTypingOutput:
        """
        Runs the sequence typing workflow using KMA.
        :param fastq_input: FASTQ input
        :param db_key: Database key
        :param db_path: Database path
        :return: None
        """
        wrapper = SequenceTypingWrapper(self._args.working_dir)
        workflow_input = SequenceTypingInput(
            fasta=ToolIOFile(self._args.fasta) if self._args.fasta else None,
            sample_name=self._sample_name,
            fastq=fastq_input,
            db_key=db_key,
            db_path=db_path,
            input_type=self._args.input_type,
        )
        wrapper.run_workflow_kma(workflow_input, self._args.threads)
        return wrapper.output

    def __export_output(self, output: SequenceTypingOutput, report: HtmlReport) -> None:
        """
        Exports the output of the workflow.
        :param output: Output
        :return: None
        """
        self._helper.logs['typing'] = output.log_file
        self._helper.informs.extend(output.informs)
        self._helper.export_output_and_commands_section(report, output.report_section)
        if self._args.output_tsv is not None:
            if output.tsv is None:
                raise ValueError("Cannot create TSV output for mixed schemes (DNA + peptide loci)")
            shutil.copyfile(output.tsv, self._args.output_tsv)


def run(args: Sequence[str] | None = None) -> None:
    """
    Entry point for the common interface.
    :param args: Command line arguments
    :return: None
    """
    script = MainSequenceTyping(args)
    script.run()


if __name__ == '__main__':
    Camel.get_instance()
    run()
