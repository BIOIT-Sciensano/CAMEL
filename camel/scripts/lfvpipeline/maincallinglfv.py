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
from camel.app.scriptutils.basepipe import basepipeutils
from camel.app.scriptutils.basepipe.basepipeutils import dict_merge
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.scripts.lfvpipeline.snakefile import variant_calling_lofreq
from camel.scripts.lfvpipeline.snakefile.variant_calling_lofreq import SNAKEFILE


@dataclasses.dataclass(frozen=True)
class Options(model.BaseOptions):
    """
    Defines the custom options for the script.
    """
    reference: Path = dataclasses.field(metadata={'help': 'Reference FASTA file'})
    output_vcf: Path = dataclasses.field(default='output.vcf', metadata={'help': 'Output VCF file'})
    gff: Path = dataclasses.field(default=None, metadata={'help': 'GFF file'})
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
    min_af: float = dataclasses.field(default=0.0, metadata={'help': 'Minimum allele frequency'})
    report_include_bam: bool = dataclasses.field(
        default=False, metadata={"help": "Include the BAM file in the output report"}
    )


class MainCalling(BaseScript[ScriptInput, ScriptOutput, Options]):
    """
    Class to perform variant calling from a BAM file.
    """

    VARIANT_CALLING_OPTS = ['call_indels', 'only_indels', 'min_af', 'report_include_bam']

    def __init__(self, in_: ScriptInput, out_: ScriptOutput, opts: Options) -> None:
        """
        Initializes the main script.
        :param in_: Input files
        :param out_: Output VCF file
        :param opts: Additional options
        :return: None
        """
        super().__init__(
            name='Low-frequency variant detection pipeline',
            version='1.0',
            title=None,
            script_in=in_,
            script_out=out_,
            script_opts=opts
        )

    def _symlink_input_data(self) -> None:
        """
        Symlinks input data (fastq and fasta)
        :return: None
        """
        dir_symlinks = self._script_opts.working_dir / 'input'
        dir_symlinks.mkdir(parents=True, exist_ok=True)

        # Script input
        self._script_in = self._script_in.symlink(dir_symlinks)

        # Reference genome
        path_ref = dir_symlinks / f'{self._script_opts.name}'
        to_replace = {'reference': path_ref}
        if not path_ref.is_symlink():
            path_ref.symlink_to(self._script_opts.reference.resolve().absolute())

        # GFF file (if provided)
        if self._script_opts.gff is not None:
            to_replace['gff'] = dir_symlinks / self._script_opts.gff.name
            gff_link = to_replace['gff']
            if not gff_link.is_symlink():
                gff_link.symlink_to(self._script_opts.gff.resolve().absolute())

        self._script_opts = dataclasses.replace(self._script_opts, **to_replace)

    def _execute(self) -> None:
        """
        Runs the variant calling Snakefile to call the variants.
        :return: None
        """
        # Symlink the fastq input
        self._symlink_input_data()
        basepipeutils.prepare_galaxy_output(self._script_out.dir, self._script_out.html)

        # Create the config file
        config_data = self.__create_snakemake_config_data()
        config_file = snakepipelineutils.generate_config_file(config_data, self._script_opts.working_dir)
        self._script_out.dir.mkdir(parents=True, exist_ok=True)

        # Run Snakemake to generate the output file
        snakepipelineutils.run_snakemake(
            snakefile=SNAKEFILE,
            config_path=config_file,
            targets=[Path(self._script_out.html), Path(variant_calling_lofreq.OUTPUT_UNFILTERED_VCF)],
            working_dir=self._script_opts.working_dir,
            threads=self._script_opts.threads)

        # Copy output
        logger.info("Collecting Snakemake output file")
        output_vcf_path = snakemakeutils.load_object(
            self._script_opts.working_dir / variant_calling_lofreq.OUTPUT_UNFILTERED_VCF)[0].path
        shutil.copyfile(output_vcf_path, self._script_opts.output_vcf)

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
                'title': self.title
            },
            'reference': {
                'name': self._script_opts.reference.name,
                'fasta': str(self._script_opts.reference),
                'gff': str(self._script_opts.gff) if self._script_opts.gff is not None else None
            },
            'variant_calling': {},
            'variant_filtering': {},
            'working_dir': str(self._script_opts.working_dir),
        }
        variant_calling_config = {
            k: getattr(self._script_opts, k)
            for k in MainCalling.VARIANT_CALLING_OPTS
            if getattr(self._script_opts, k, None) is not None
        }
        if variant_calling_config['only_indels']:
            variant_calling_config['call_indels'] = True

        config_data['variant_calling'].update(variant_calling_config)
        if self._script_opts.__getattribute__('gff') is not None:
            config_data['variant_calling']['csq'] = True
        return config_data


@click.command(name='lfv_pipeline', short_help='Detection of low frequency variants.')
@basescriptutils.add_input_opts()
@basescriptutils.add_output_opts
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Entry point for the common interface.
    :param kwargs: Command line arguments
    :return: None
    """
    script_input = basescriptutils.parse_script_input(kwargs)
    script_out = basescriptutils.parse_script_output(kwargs)

    script = MainCalling(
        in_=script_input,
        out_=script_out,
        opts=Options(**cliutils.from_kwargs(Options, kwargs))
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
