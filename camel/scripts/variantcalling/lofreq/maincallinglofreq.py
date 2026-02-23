#!/usr/bin/env python

import dataclasses
import shutil
from pathlib import Path
from typing import Any

import click

from camel.app.cli import cliutils
from camel.app.core.snakemake import snakepipelineutils, snakemakeutils
from camel.app.loggers import logger, initialize_logging
from camel.app.scriptutils import model
from camel.app.scriptutils.basepipe.basepipeutils import dict_merge
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.scripts.variantcalling.lofreq.snakefile import variant_calling_lofreq
from camel.scripts.variantcalling.lofreq.snakefile.variant_calling_lofreq import SNAKEFILE


@dataclasses.dataclass(frozen=True)
class Output(model.BaseOutput):
    """
    Defines the script output.
    """
    output: Path = dataclasses.field(metadata={'help': 'Output VCF file'})
    html: Path = dataclasses.field(metadata={'help': 'Output HTML report'})
    dir: Path = dataclasses.field(metadata={'help': 'Output directory'})


@dataclasses.dataclass(frozen=True)
class Options(model.BaseOptions):
    """
    Defines the custom options for the script.
    """
    reference: Path = dataclasses.field(metadata={'help': 'Reference FASTA file'})
    reference_name: str | None = dataclasses.field(default=None, metadata={'help': 'Reference genome name'})

    @property
    def name(self) -> str:
        """
        Returns the input name.
        :return: Input name
        """
        name = self.reference_name if self.reference_name is not None else self.reference.name
        return name

    working_dir: Path = dataclasses.field(default=Path.cwd(), metadata={'help': 'Working directory'})
    call_indels: bool = dataclasses.field(default=False, metadata={'help': 'Call indels'})
    only_indels: bool = dataclasses.field(default=False, metadata={'help': 'Call only indels'})
    threads: int = dataclasses.field(default=8, metadata={'help': 'Number of threads'})


class MainCalling(BaseScript[ScriptInput, Output, Options]):
    """
    Class to perform variant calling from a BAM file.
    """

    def __init__(self, in_: ScriptInput, out_: Output, opts: Options) -> None:
        """
        Initializes the main script.
        :param in_: Input files
        :param out_: Output VCF file
        :param opts: Additional options
        :return: None
        """
        super().__init__(
            name='Variant calling (LoFreq)',
            version='1.0',
            script_in=in_,
            script_out=out_,
            script_opts=opts
        )

    def _symlink_input_data(self):
        """
        Symlinks input data (fastq and fasta)
        :return: None
        """
        input_symlink_dir = self._script_opts.working_dir / 'input'
        input_symlink_dir.mkdir(parents=True, exist_ok=True)

        links = []

        links.extend([('fastq_pe', self._script_in.fastq_pe[0],
                       self._script_in.fastq_pe_names[0] if self._script_in.fastq_pe_names else
                       self._script_in.fastq_pe[0].name)])
        links.extend([('fastq_pe', self._script_in.fastq_pe[1],
                       self._script_in.fastq_pe_names[1] if self._script_in.fastq_pe_names else
                       self._script_in.fastq_pe[1].name)])
        links.extend([('fasta', self._script_opts.reference,
                       self._script_opts.reference_name if self._script_opts.reference_name else
                       self._script_opts.reference.name)])

        to_replace = {}
        for key, origin, target in links:
            path_out = input_symlink_dir / target
            if path_out.is_symlink():
                logger.debug(f'Symlink already exists: {path_out}')
            else:
                ((input_symlink_dir / target).absolute()).symlink_to(origin.resolve().absolute())

            # Save the updated path
            if key not in to_replace:
                to_replace[key] = path_out
            else:
                to_replace[key] = [to_replace[key], path_out]

        # Update the script input
        script_in_updated = dataclasses.replace(self._script_in, **to_replace)
        self._script_in = script_in_updated

    def _execute(self) -> None:
        """
        Runs the variant calling Snakefile to call the variants.
        :return: None
        """
        # Symlink the fastq input
        self._symlink_input_data()

        # Create the config file
        config_data = self.__create_snakemake_config_data()
        config_file = snakepipelineutils.generate_config_file(config_data, self._script_opts.working_dir)

        # Run Snakemake to generate the output file
        snakepipelineutils.run_snakemake(
            snakefile=SNAKEFILE,
            config_path=config_file,
            targets=[Path(variant_calling_lofreq.OUTPUT_REPORT)],
            working_dir=self._script_opts.working_dir,
            threads=self._script_opts.threads)

        # Copy output
        logger.info("Collecting Snakemake output file")
        output_vcf_path = snakemakeutils.load_object(self._script_opts.working_dir /
                                                     variant_calling_lofreq.OUTPUT_UNFILTERED_VCF)[0].path
        shutil.copyfile(output_vcf_path, self._script_out.output)
        if not self._script_out.html.parent.exists():
            self._script_out.html.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(self._script_opts.working_dir / variant_calling_lofreq.OUTPUT_REPORT, self._script_out.html)

    def __create_snakemake_config_data(self) -> dict:
        """
        Creates a Snakemake configuration file.
        :return: Config file data
        """
        input_dict = self._script_in.to_dict()
        dict_merge(input_dict, {'name': self._script_in.name,
                                'input_str': self._script_in.input_str
                                })
        config_data: dict[str, Any] = {
            'input': input_dict,
            'output': {
                'dir': str(self._script_out.dir),
                'html': str(self._script_out.html)
            },
            'read_trimming': {
                'method': 'fastp'
            },
            'script_info': {
                'name': self.name,
                'version': self.version,
            },
            'reference': {
                'name': self._script_opts.reference.name,
                'fasta': str(self._script_opts.reference),
            },
            'variant_calling': {},
            'variant_filtering': {},
            'working_dir': str(self._script_opts.working_dir),
        }
        for k in ['call_indels', 'only_indels']:
            if self._script_opts.__getattribute__(k) is not None:
                config_data['variant_calling'][k] = self._script_opts.__getattribute__(k)
        return config_data


@click.command(name='variant_calling_lofreq', short_help='Call variants from FASTQ files using LoFreq')
@basescriptutils.add_input_opts()
@cliutils.add_click_options_from_dataclass(Output)
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Entry point for the common interface.
    :param kwargs: Command line arguments
    :return: None
    """
    script_input = basescriptutils.parse_script_input(kwargs)

    script = MainCalling(
        in_=script_input,
        out_=Output(**cliutils.from_kwargs(Output, kwargs)),
        opts=Options(**cliutils.from_kwargs(Options, kwargs))
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
