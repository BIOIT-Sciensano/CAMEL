#!/usr/bin/env python
import dataclasses
from pathlib import Path

import click
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.cli import cliutils
from camel.app.loggers import initialize_logging
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.model import BaseInput, BaseOptions, BaseOutput
from camel.app.tools.freebayes.freebayes import Freebayes


@dataclasses.dataclass(frozen=True)
class Input(BaseInput):
    """
    Input for the FreeBayes script.
    """
    bam: Path = dataclasses.field(metadata={'help': 'Input BAM file'})
    reference: Path = dataclasses.field(metadata={'help': 'Reference FASTA file'})

@dataclasses.dataclass(frozen=True)
class Output(BaseOutput):
    """
    Output for the FreeBayes script.
    """
    output: Path = dataclasses.field(metadata={'help': 'Output VCF file'})

@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Options for the FreeBayes script.
    """
    working_dir: Path = dataclasses.field(default=Path.cwd(), metadata={'help': 'Working directory'})
    min_base_quality: int = dataclasses.field(default=0, metadata={'help': 'Minimum base quality'})
    min_coverage: int = dataclasses.field(default=0, metadata={'help': 'Minimum coverage'})
    min_mapping_quality: int = dataclasses.field(default=1, metadata={'help': 'Minimum mapping quality'})
    min_supporting_allele_qsum: int = dataclasses.field(default=0, metadata={'help': 'Minimum supporting allele QSum'})
    ploidy: int = dataclasses.field(default=1, metadata={'help': 'Ploidy'})
    report_monomorphic: bool = dataclasses.field(default=False, metadata={'help': 'Report monomorphic sites'})
    standard_filters: bool = dataclasses.field(default=False, metadata={'help': 'Use standard filters'})


class MainFreebayesCalling(BaseScript[Input, Output, Options]):
    """
    Main script for the FreeBayes tool.
    """

    def __init__(self, in_: Input, out: Output, opts: Options) -> None:
        """
        Initializes the main script.
        :param in_: Script input
        :param out: Script output
        :param opts: Options
        :return: None
        """
        super().__init__(
            name='Freebayes',
            version='1.0.0',
            script_in=in_,
            script_out=out,
            script_opts=opts
        )

    def _execute(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        freebayes = Freebayes()

        # Input
        if self._script_in.reference is not None:
            freebayes.add_input_files({'FASTA': [ToolIOFile(self._script_in.reference)]})
        if self._script_in.bam is not None:
            freebayes.add_input_files({'BAM': [ToolIOFile(self._script_in.bam)]})

        # Options
        if self._script_opts.standard_filters:
            freebayes.update_parameters(standard_filters=True)
        else:
            if self._script_opts.min_base_quality:
                freebayes.update_parameters(min_base_quality=self._script_opts.min_base_quality)
            if self._script_opts.min_mapping_quality != 1:
                freebayes.update_parameters(min_mapping_quality=self._script_opts.min_mapping_quality)
            if self._script_opts.min_supporting_allele_qsum:
                freebayes.update_parameters(min_supporting_allele_qsum=self._script_opts.min_supporting_allele_qsum)

        if self._script_opts.report_monomorphic:
            freebayes.update_parameters(report_monomorphic=True)
        if self._script_opts.min_coverage:
            freebayes.update_parameters(min_coverage=self._script_opts.min_coverage)
        freebayes.update_parameters(ploidy=self._script_opts.ploidy, vcf=str(self._script_out.output))

        # Run the tool
        freebayes.run(self._script_opts.working_dir)

@click.command(name='freebayes', short_help='Variant calling using freebayes')
@cliutils.add_click_options_from_dataclass(Input)
@cliutils.add_click_options_from_dataclass(Output)
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Variant calling using freebayes.
    """
    script = MainFreebayesCalling(
        in_=Input(**cliutils.from_kwargs(Input, kwargs)),
        out=Output(**cliutils.from_kwargs(Output, kwargs)),
        opts=Options(**cliutils.from_kwargs(Options, kwargs))
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
