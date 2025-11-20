#!/usr/bin/env python
import dataclasses
import shutil
from pathlib import Path
from typing import Any

import click

from camel.app.cli import cliutils
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.snakemake import snakemakeutils
from camel.app.core.snakemake import snakepipelineutils
from camel.app.loggers import logger, initialize_logging
from camel.app.scriptutils import model
from camel.app.scriptutils.basescript.bamwithrefinput import BAMWithRefInput
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.snakefiles import variant_calling


@dataclasses.dataclass(frozen=True)
class Output(model.BaseOutput):
    """
    Defines the script output.
    """
    output: Path = dataclasses.field(metadata={'help': 'Output BAM file'})
    output_consensus: Path | None = dataclasses.field(metadata={'help': 'Output consensus FASTA file'})

@dataclasses.dataclass(frozen=True)
class Options(model.BaseOptions):
    """
    Defines the custom options for the script.
    """
    working_dir: Path = dataclasses.field(default=Path.cwd(), metadata={'help': 'Working directory'})
    input_type: str = dataclasses.field(default=model.InputType.ILLUMINA.value, metadata={
        'help': 'Input type',
        'choices': [model.InputType.ILLUMINA.value, model.InputType.ONT.value]
    })
    ploidy: str = dataclasses.field(default='1', metadata={
        'help': 'Ploidy', 'choices': ['GRCh37', 'GRCh38', 'X', 'Y', '1']})
    calling_method: str = dataclasses.field(default='consensus', metadata={
        'help': 'Calling method', 'choices': ['consensus', 'multiallelic']})
    skip_variants: str | None = dataclasses.field(default=None, metadata={
        'help': 'Skip variants', 'choices': ['snps', 'indels']})
    mutation_rate: float | None = dataclasses.field(default=None, metadata={'help': 'Mutation rate'})
    minimal_mq: int | None = dataclasses.field(default=None, metadata={'help': 'Minimal mapping quality'})
    minimal_bq: int | None = dataclasses.field(default=None, metadata={'help': 'Minimal base quality'})
    output_all_sites: bool = dataclasses.field(default=False, metadata={'help': 'Output all sites'})
    count_orphans: bool = dataclasses.field(default=False, metadata={'help': 'Count orphans'})
    disable_baq: bool = dataclasses.field(default=False, metadata={'help': 'Disable BAQ'})
    threads: int = dataclasses.field(default=1, metadata={'help': 'Number of threads'})

class MainCalling(BaseScript[BAMWithRefInput, Output, Options]):
    """
    Class to perform variant calling from a BAM file.
    """

    def __init__(self, in_: BAMWithRefInput, out_: Output, opts: Options) -> None:
        """
        Initializes the main script.
        """
        super().__init__(
            name='Variant calling (SAMtools)',
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
        self._script_in: BAMWithRefInput = self._script_in.create_symlinks(self._script_opts.working_dir / 'input')

        # Create the config file
        config_data = self.__create_snakemake_config_data()
        config_file = snakepipelineutils.generate_config_file(config_data, self._script_opts.working_dir)

        # Copy input the BAM file to the right location
        target_dir = self._script_opts.working_dir / 'variant_calling' / 'read_mapping' / self._script_opts.input_type
        if not target_dir.exists():
            target_dir.mkdir(parents=True)
        snakemakeutils.dump_object([ToolIOFile(self._script_in.bam)], target_dir / 'bam.io')

        # Run Snakemake to generate the output file
        snakepipelineutils.run_snakemake(
            snakefile=variant_calling.SNAKEFILE,
            config_path=config_file,
            targets=[Path(variant_calling.OUTPUT_UNFILTERED_VCF)],
            working_dir=self._script_opts.working_dir,
            threads=self._script_opts.threads)

        # Generate consensus sequence
        if self._script_out.output_consensus:
            self.__generate_consensus_sequence(self._script_out.output_consensus, config_data)

        # Copy output
        logger.info("Collecting Snakemake output file")
        output_vcf_path = snakemakeutils.load_object(self._script_opts.working_dir / variant_calling.OUTPUT_UNFILTERED_VCF)[0].path
        shutil.copyfile(output_vcf_path, self._script_out.output)

    def __create_snakemake_config_data(self) -> dict:
        """
        Creates a Snakemake configuration file.
        :return: Config file data
        """
        config_data: dict[str, Any] = {
            'input': {
                'type': self._script_opts.input_type,
                'sample_name': 'sample'
            },
            'reference': {
                'name': self._script_in.reference.name,
                'fasta': str(self._script_in.reference),
            },
            'variant_calling': {'ploidy': self._script_opts.ploidy},
            'variant_filtering': {},
            'working_dir': str(self._script_opts.working_dir),
        }
        for k in ['calling_method', 'skip_variants', 'mutation_rate', 'minimal_bq', 'minimal_mq', 'count_orphans',
                  'disable_baq']:
            if self._script_opts.__getattribute__(k) is not None:
                config_data['variant_calling'][k] = self._script_opts.__getattribute__(k)
        if self._script_opts.output_all_sites is True:
            config_data['variant_calling']['variants_only'] = False
        return config_data

    def __generate_consensus_sequence(self, output_path: Path, config_data: dict[str, Any]) -> None:
        """
        Generates the consensus sequence by applying the detected variants to the reference sequence.
        :param output_path: Output path to save the consensus sequence
        :param config_data: Snakemake config data
        :return: None
        """
        config_file = snakepipelineutils.generate_config_file(
            config_data, self._script_opts.working_dir, 'consensus.yml')
        snakepipelineutils.run_snakemake(
            variant_calling.SNAKEFILE, config_file, [Path(variant_calling.OUTPUT_CONSENSUS)],
            self._script_opts.working_dir, threads=self._script_opts.threads)
        fasta_consensus = snakemakeutils.load_object(self._script_opts.working_dir / variant_calling.OUTPUT_CONSENSUS)[0].path
        shutil.copyfile(fasta_consensus, output_path)

@click.command(name='variant_calling_samtools', short_help='Call variants from a BAM file using SAMtools')
@cliutils.add_click_options_from_dataclass(BAMWithRefInput)
@cliutils.add_click_options_from_dataclass(Output)
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Entry point for the common interface.
    :param kwargs: Command line arguments
    :return: None
    """
    script = MainCalling(
        in_=BAMWithRefInput(
            bam=kwargs['bam'],
            bam_name=kwargs['bam_name'],
            reference=kwargs['reference'],
            reference_name=kwargs['reference_name'],
        ),
        out_=Output(**cliutils.from_kwargs(Output, kwargs)),
        opts=Options(**cliutils.from_kwargs(Options, kwargs))
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
