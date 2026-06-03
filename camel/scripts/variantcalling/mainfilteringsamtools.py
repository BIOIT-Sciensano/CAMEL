#!/usr/bin/env python
import dataclasses
import shutil
from pathlib import Path
from typing import Any

import click
import yaml

from camel.app.cli import cliutils
from camel.app.loggers import initialize_logging
from camel.app.scriptutils import model
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.tools.bcftools.bcftoolsmpileup import BcftoolsMpileup
from camel.app.wrappers.variantfilteringwrapper import VariantFilteringWrapper
from camel.version import __VERSION__


@dataclasses.dataclass(frozen=True)
class Input(model.BaseInput):
    """
    Input for the variant filtering.
    """
    vcf: Path = dataclasses.field(metadata={'help': 'Input VCF file'})
    bam: Path | None = dataclasses.field(metadata={'help': 'Input BAM file'})

    def name(self) -> str:
        """
        Returns the dataset name.
        :return: Name
        """
        return self.bam.name


@dataclasses.dataclass(frozen=True)
class Output(model.BaseOutput):
    """
    Output for the variant filtering.
    """
    output_vcf: Path = dataclasses.field(metadata={'help': 'Output VCF file'})
    output_stats: Path | None = dataclasses.field(metadata={'help': 'Filtering statistics (in JSON format)'})


@dataclasses.dataclass(frozen=True)
class Options(model.BaseOptions):
    """
    Options for the variant filtering.
    """
    working_dir: Path = dataclasses.field(default=Path.cwd(), metadata={'help': 'Working directory'})
    input_type: str = dataclasses.field(
        default=model.InputType.ILLUMINA.value,
        metadata={
            'choices': [model.InputType.ILLUMINA.value, model.InputType.ONT.value],
            'help': 'Input type'}
    )
    bed: Path | None = dataclasses.field(default=None, metadata={'help': 'BED file with regions to remove'})
    keep_best: bool = dataclasses.field(default=False, metadata={'help': 'Keep best variant'})
    min_distance: int = dataclasses.field(default=10, metadata={'help': 'Minimum distance between variants'})
    min_forward_depth: int = dataclasses.field(default=1, metadata={'help': 'Minimum forward depth'})
    min_mapping_quality: int = dataclasses.field(default=30, metadata={'help': 'Minimum mapping quality'})
    min_reverse_depth: int = dataclasses.field(default=1, metadata={'help': 'Minimum reverse depth'})
    min_snp_quality: float = dataclasses.field(default=25, metadata={'help': 'Minimum SNP quality'})
    min_total_depth: int = dataclasses.field(default=10, metadata={'help': 'Minimum total depth'})
    min_zscore: float = dataclasses.field(default=1.96, metadata={'help': 'Minimum Z-score'})
    soft_filter: bool = dataclasses.field(default=False, metadata={'help': 'Soft filter'})
    y_mult: float = dataclasses.field(default=10, metadata={'help': 'Y multiplier'})

class MainFiltering(BaseScript[Input, Output, Options]):
    """
    Class to run the samtools variant filtering using CAMEL.
    """

    def __init__(self, script_in: Input, script_out: Output, opts: Options) -> None:
        """
        Initializes the main script.
        :param script_in: Script input
        :param script_out: Script output
        :param opts: Script options
        :return: None
        """
        tool_version = BcftoolsMpileup().version
        super().__init__(
            name='Variant filtering (SAMtools)',
            version=f'{tool_version}+CAMEL_{__VERSION__}',
            script_in=script_in,
            script_out=script_out,
            script_opts=opts
        )

    def _execute(self) -> None:
        """
        Filters the input VCF file.
        :return: None
        """
        wrapper = VariantFilteringWrapper(self._script_opts.working_dir)
        wrapper.run(
            sample_name=self._script_in.vcf.stem, vcf_file=self._script_in.vcf, bam_file=self._script_in.bam,
            filtering_options=self.__get_filtering_options(), input_type=self._script_opts.input_type)
        shutil.copyfile(wrapper.output.vcf_filtered.path, self._script_out.output_vcf)
        if self._script_out.output_stats is not None:
            with open(self._script_out.output_stats, 'w') as handle:
                yaml.dump(wrapper.output.stats, handle)

    def __get_filtering_options(self) -> dict[str, Any]:
        """
        Returns the dictionary with filtering options.
        :return: Filtering options
        """
        filtering_opts = {
            'soft_filter': self._script_opts.soft_filter,
            'input': {'type': self._script_opts.input_type},
            'depth': {
                'min_total_depth': self._script_opts.min_total_depth,
                'min_fwd_depth': self._script_opts.min_forward_depth,
                'min_rev_depth': self._script_opts.min_reverse_depth},
            'snp_quality': {
                'min_snp_quality': self._script_opts.min_snp_quality},
            'mapping_quality': {
                'min_mapping_quality': self._script_opts.min_mapping_quality},
            'distance': {
                'min_distance': self._script_opts.min_distance,
                'keep_best': self._script_opts.keep_best},
            'zscore': {
                'min_zscore': self._script_opts.min_zscore,
                'y_multiplier': self._script_opts.y_mult},
        }
        if self._script_opts.bed is not None:
            filtering_opts['region'] = {'bed_file': str(self._script_opts.bed)}
        return filtering_opts


@click.command(name='variant_filtering_samtools', short_help='Filter variants from a VCF file called using SAMtools')
@cliutils.add_click_options_from_dataclass(Input)
@cliutils.add_click_options_from_dataclass(Output)
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Filter variants from a VCF file called using SAMtools.
    """
    script = MainFiltering(
        script_in=Input(**cliutils.from_kwargs(Input, kwargs)),
        script_out=Output(**cliutils.from_kwargs(Output, kwargs)),
        opts=Options(**cliutils.from_kwargs(Options, kwargs))
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
