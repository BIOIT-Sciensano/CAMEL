#!/usr/bin/env python
import dataclasses
import shutil
from pathlib import Path

import click

from camel.app.cli import cliutils
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.reports import reportutils
from camel.app.core.reports.htmlreport import HtmlReport
from camel.app.loggers import initialize_logging
from camel.app.scriptutils import inputhelper, model
from camel.app.scriptutils.basepipe.fastqinput import FastqInput
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.inputhelper import helper_by_input_type
from camel.app.scriptutils.inputhelper.inputhelperbase import InputHelperBase
from camel.app.scriptutils.model import BaseOptions, BaseOutput
from camel.app.toolkits.sequencetyping.sequencetypingutils import SequenceTypingUtils
from camel.app.wrappers.sequencetypingwrapper import (
    SequenceTypingInput,
    SequenceTypingOutput,
    SequenceTypingWrapper,
)


@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Custom options for the sequence typing script.
    """
    db: Path = dataclasses.field(default=None, metadata={'help': 'Path to sequence typing database'})
    detection_method: str = dataclasses.field(default='blast', metadata={
        'help': 'Sequence typing method',
        'choices': ['blast', 'kma', 'mist'],
        'show_default': True})
    working_dir: Path = dataclasses.field(default=Path('work'), metadata={
        'help': 'Working directory',
        'show_default': True})
    threads: int | None = dataclasses.field(default=1, metadata={'help': 'Number of threads', 'show_default': True})
    blastn_task: str = dataclasses.field(default='megablast', metadata={
        'help': 'BLASTn task', "show_default": True, 'choices': ['blastn', 'megablast']})

@dataclasses.dataclass(frozen=True)
class Output(BaseOutput):
    """
    Output for the sequence typing script.
    """
    output_html: Path = dataclasses.field(metadata={'help': 'Output HTML file'})
    output_dir: Path | None = dataclasses.field(metadata={'help': 'Output directory'})
    output_tsv: Path | None = dataclasses.field(default=None, metadata={'help': 'Output TSV file'})
    output_fasta: Path | None = dataclasses.field(default=None, metadata={
        'help': 'Output FASTA file (if reads are assembled)'})


class MainSequenceTyping(BaseScript[ScriptInput, Output, Options]):
    """
    Main script for the sequence typing workflow.
    """

    def __init__(self, script_in: ScriptInput, helper: InputHelperBase, script_out: Output, script_opts: Options) -> None:
        """
        Initializes the main script.
        :param script_in: Script input
        :param helper: Input helper
        :param script_out: Script output
        :param script_opts: Script options
        """
        super().__init__(
            name='Sequence typing',
            version='1.0',
            script_in=script_in,
            script_out=script_out,
            script_opts=script_opts
        )
        self._helper: InputHelperBase = helper

    def _execute(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        # Initialize report
        report = reportutils.init_report(
            path_out=self._script_out.output_html,
            dir_out=self._script_out.output_dir,
            key='Sequence typing',
            title=f'Sequence typing {self._script_opts.detection_method}'
        )
        report.add_html_object(reportutils.create_overview_section(
            version=self._version,
            dataset_name=self._script_in.name,
            input_file_str=self._script_in.input_str,
            extra_data=[('Detection method', self._script_opts.detection_method)]
        ))
        report.save()

        # Run script with wrapper
        db_data = SequenceTypingUtils.parse_scheme_metadata(self._script_opts.db)
        if self._script_opts.detection_method in ('blast', 'mist'):
            path_fasta = self._helper.prepare_fasta_input(self._script_in, report)
            # Save assembly if specified
            if self._script_out.output_fasta is not None:
                shutil.copyfile(str(path_fasta), self._script_out.output_fasta)
            if self._script_opts.detection_method == 'blast':
                output = self.__run_sequence_typing_blast(path_fasta, db_data['name'], self._script_opts.db)
            else:
                output = self.__run_sequence_typing_mist(path_fasta, db_data['name'], self._script_opts.db)
        elif self._script_opts.detection_method == 'kma':
            fastq_input = self._helper.prepare_fastq_input(self._script_in, report)
            output = self.__run_sequence_typing_kma(fastq_input, db_data['name'], self._script_opts.db)
        else:
            raise ValueError(f"Invalid detection method: {self._script_opts.detection_method}")
        self.__export_output(output, report)

    def __run_sequence_typing_blast(self, fasta_file: Path, db_key: str, db_path: Path) -> SequenceTypingOutput:
        """
        Runs the sequence typing workflow using BLAST.
        :param fasta_file: Input FASTA file
        :param db_key: Database key
        :param db_path: Database directory path
        :return: None
        """
        wrapper = SequenceTypingWrapper(self._script_opts.working_dir)
        workflow_input = SequenceTypingInput(
            sample_name=self._script_in.name,
            fasta=ToolIOFile(fasta_file),
            input_type='fasta',
            db_path=db_path,
            db_key=db_key
        )
        wrapper.run_workflow_blast(workflow_input, self._script_opts.blastn_task, self._script_opts.threads)
        return wrapper.output

    def __run_sequence_typing_mist(self, fasta_file: Path, db_key: str, db_path: Path) -> SequenceTypingOutput:
        """
        Runs the sequence typing workflow using MiST.
        :param fasta_file: Input FASTA file
        :param db_key: Database key
        :param db_path: Database directory path
        :return: None
        """
        wrapper = SequenceTypingWrapper(self._script_opts.working_dir)
        workflow_input = SequenceTypingInput(
            sample_name=self._script_in.name,
            fasta=ToolIOFile(fasta_file),
            input_type='fasta',
            db_path=db_path,
            db_key=db_key
        )
        wrapper.run_workflow_mist(workflow_input, threads=self._script_opts.threads)
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
        wrapper = SequenceTypingWrapper(self._script_opts.working_dir)
        workflow_input = SequenceTypingInput(
            fasta=ToolIOFile(self._script_in.fasta) if self._script_in.fasta else None,
            sample_name=self._script_in.name,
            fastq=fastq_input,
            db_key=db_key,
            db_path=db_path,
            input_type=self._script_in.type_.value,
        )
        wrapper.run_workflow_kma(workflow_input, self._script_opts.threads)
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
        if self._script_out.output_tsv is not None:
            if output.tsv is None:
                raise ValueError("Cannot create TSV output for mixed schemes (DNA + peptide loci)")
            shutil.copyfile(output.tsv, self._script_out.output_tsv)


@click.command(name='sequence_typing', short_help='Sequence typing (MLST, cgMLST, etc)')
@basescriptutils.add_input_opts(supported=[model.InputType.FASTA, model.InputType.ILLUMINA, model.InputType.ONT])
@cliutils.add_click_options_from_dataclass(Output)
@cliutils.add_click_options_from_dataclass(Options)
@inputhelper.add_helper_opts
def main(**kwargs) -> None:
    """
    Runs the script.
    :return: None
    """
    # Parse the script input
    script_input = basescriptutils.parse_script_input(kwargs)
    script_opts = Options(**cliutils.from_kwargs(Options, kwargs))

    # Initialize the helper class to prepare the input
    helper = helper_by_input_type[script_input.type_](dir_=script_opts.working_dir, name=script_input.name)
    helper.set_opts(*helper.opts_from_cli(kwargs))

    # Run the main script
    script = MainSequenceTyping(
        script_in=basescriptutils.parse_script_input(kwargs),
        script_out=Output(**cliutils.from_kwargs(Output, kwargs)),
        script_opts=script_opts,
        helper=helper
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
