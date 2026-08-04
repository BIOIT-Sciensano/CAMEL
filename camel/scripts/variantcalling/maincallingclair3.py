#!/usr/bin/env python
import dataclasses
import shutil
from pathlib import Path

import click

from camel.app.cli import cliutils
from camelcore.app.io.tooliofile import ToolIOFile
from camel.app.core.snakemake import snakemakeutils
from camel.app.core.snakemake import snakepipelineutils
from camel.app.loggers import logger, initialize_logging
from camel.app.scriptutils import model
from camel.app.scriptutils.basescript.bamwithrefinput import BAMWithRefInput
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.snakefiles import variant_calling_clair3


@dataclasses.dataclass(frozen=True)
class Output(model.BaseOutput):
    """
    Defines the script output.
    """
    output: Path = dataclasses.field(metadata={'help': 'Output BAM file'})

@dataclasses.dataclass(frozen=True)
class Options(model.BaseOptions):
    """
    Defines the custom options for the script.
    """
    model_path: Path = dataclasses.field(metadata={'help': 'Path to the model directory'})
    working_dir: Path = dataclasses.field(default=Path.cwd(), metadata={'help': 'Working directory'})
    platform: str = dataclasses.field(default='ilmn', metadata={'choices': ['ont', 'hifi', 'ilmn']})
    haploid_precise: bool = dataclasses.field(default=False, metadata={'help': 'Use haploid-precise mode'})
    no_phasing: bool = dataclasses.field(default=False, metadata={'help': 'Disable phasing'})
    include_ctgs: bool = dataclasses.field(default=False, metadata={'help': 'Include contigs'})
    long_indel: bool = dataclasses.field(default=False, metadata={'help': 'Use long indel mode'})
    threads: int = dataclasses.field(default=8, metadata={'help': 'Number of threads'})

class MainCallingClair3(BaseScript[BAMWithRefInput, Output, Options]):
    """
    Class to run Clair3 variant calling using CAMEL.
    """

    def __init__(self, in_: BAMWithRefInput, out_: Output, opts: Options) -> None:
        """
        Initializes the main script.
        """
        super().__init__(
            name='Variant calling (Clair3)',
            version='1.0',
            script_in=in_,
            script_out=out_,
            script_opts=opts
        )

    def _execute(self) -> None:
        """
        Runs the variant calling Snakefile to call the variants.
        :return: None
        """
        # Create the configuration file
        config_data = self.__create_snakemake_config_data()
        config_file = snakepipelineutils.generate_config_file(config_data, self._script_opts.working_dir)

        # Copy input BAM file to the right location
        target_dir = self._script_opts.working_dir / 'variant_calling' / 'read_mapping'
        if not target_dir.exists():
            target_dir.mkdir(parents=True)
        snakemakeutils.dump_object([ToolIOFile(self._script_in.bam)], target_dir / 'bam.io')

        # Run Snakemake to generate output file
        path_vcf = variant_calling_clair3.OUTPUT_UNFILTERED_VCF
        snakepipelineutils.run_snakemake(
            variant_calling_clair3.SNAKEFILE, config_file, [path_vcf], self._script_opts.working_dir,
            self._script_opts.threads)

        # Copy output
        logger.info("Collecting Snakemake output file")
        output_vcf_path = snakemakeutils.load_object(self._script_opts.working_dir / path_vcf)[0].path
        shutil.copyfile(output_vcf_path, self._script_out.output)

    def __create_snakemake_config_data(self) -> dict:
        """
        Creates a Snakemake configuration file.
        :return: Config file data
        """
        config_data = {
            'input': {'sample_name': self._script_in.name},
            'working_dir': str(self._script_opts.working_dir),
            'variant_calling': {
                'platform': self._script_opts.platform,
                'reference': {
                    'name': self._script_in.reference_name if self._script_in.reference_name else self._script_in.reference.name,
                    'path': str(self._script_in.reference)},
                'bam': str(self._script_in.bam),
                'clair3': {
                    'model_path': str(self._script_opts.model_path),
                    'haploid_precise': self._script_opts.haploid_precise,
                    'no_phasing': self._script_opts.no_phasing,
                    'include_ctgs': self._script_opts.include_ctgs,
                    'long_indel': self._script_opts.long_indel,
                    'platform': self._script_opts.platform
                }
            },
        }
        return config_data


@click.command(name='variant_calling_clair', short_help='Call variants from a BAM file using Clair3')
@cliutils.add_click_options_from_dataclass(BAMWithRefInput)
@cliutils.add_click_options_from_dataclass(Output)
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Entry point for the common interface.
    :param kwargs: Command line arguments
    :return: None
    """
    script = MainCallingClair3(
        in_=BAMWithRefInput(**cliutils.from_kwargs(BAMWithRefInput, kwargs)),
        out_=Output(**cliutils.from_kwargs(Output, kwargs)),
        opts=Options(**cliutils.from_kwargs(Options, kwargs))
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
