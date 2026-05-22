#!/usr/bin/env python
import dataclasses
import json
import re
import shutil
from importlib.resources import files
from pathlib import Path

import click
import yaml
from Bio import SeqIO
from camelcore.app.utils import fastautils

from camel.app.cli import cliutils
from camel.app.config import config
from camel.app.core.snakemake import snakepipelineutils, snakemakeutils
from camel.app.loggers import logger, initialize_logging
from camel.app.scriptutils import model
from camel.app.scriptutils.basepipe import basepipeutils
from camel.app.scriptutils.basepipe.basepipe import BasePipe
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoptions import ScriptOptions
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.scripts.viralconsensuspipeline import SNAKEFILE_MAIN
from camel.scripts.viralconsensuspipeline.snakefile import iterativemapping


DB_ROOT = Path(config.dir_db, 'pipelines', 'viral_consensus', 'version_1.1')

@dataclasses.dataclass(frozen=True)
class Options(model.BaseOptions):
    """
    Pipeline-specific options.
    """
    # Species
    species: str
    species_name: str | None = None

    # Reference genome or database
    fasta_ref: Path | None = dataclasses.field(default=None, metadata={'help': 'Reference genome FASTA file'})
    fasta_ref_name: str | None = None
    ref_genome_db: Path | None = dataclasses.field(default=None, metadata={
        'help': 'Path to reference genome database'})

    # Primer removal
    fasta_primers: Path | None = dataclasses.field(default=None, metadata={
        'help': 'Path to FASTA file with primer sequences'})
    fasta_primers_name: str | None = None

    # Downsampling & gaps
    cov_max_segment: int = dataclasses.field(default=10_000, metadata={'help': 'Maximum segment coverage'})
    gap_depth_cutoff: int = dataclasses.field(default=5, metadata={
        'help': 'Positions with a depth smaller than this value are flagged as missing / gaps'})
    gap_len_cutoff: int = dataclasses.field(default=10, metadata={
        'help': 'Minimum length to mark a region as a gap'})

    # Variant filtering & iterative mapping
    max_iter: int = dataclasses.field(default=6, metadata={
        'help': 'Maximum number of iterations for iterative mapping'})
    variant_min_af: float = dataclasses.field(default=0.5, metadata={'help': 'Minimum allele frequency'})
    variant_min_dp: float = dataclasses.field(default=5, metadata={'help': 'Minimum depth at variant position'})
    variant_min_qual: int = dataclasses.field(default=10, metadata={'help': 'Minimum variant quality'})
    variant_min_mq: int = dataclasses.field(default=30, metadata={'help': 'Minimum mapping quality'})
    clair3_model: Path | None = dataclasses.field(default=None, metadata={'help': 'Clair3 variant calling model'})

    # Other options
    analyses: list[str] = dataclasses.field(default_factory=list)


class MainViralConsensusPipeline(BasePipe):
    """
    Main script for the viral consensus pipeline.
    """

    def __init__(
        self,
        in_: ScriptInput,
        out: ScriptOutput,
        opts: ScriptOptions,
        opts_custom: Options
    ) -> None:
        """
        Initializes the main class.
        :param in_: Script input
        :param out: Script output
        :param opts: General pipeline options
        :param opts_custom: Pipeline-specific options
        :return: None
        """
        super().__init__(
            name='Viral consensus pipeline',
            version='1.1',
            script_in=in_,
            script_out=out,
            opts=opts,
            snakefile=SNAKEFILE_MAIN
        )
        self._opts_custom = opts_custom

    def _validate_opts(self) -> None:
        """
        Checks if the provided options are valid.
        :return: None
        """
        # Species name
        if (self._opts_custom.species == 'other') and (self._opts_custom.species_name is None):
            raise ValueError("Species name needs to be specified when species is 'other'")

        # Reference genome
        if self._opts_custom.fasta_ref is not None:
            MainViralConsensusPipeline._validate_ref_genome_file(self._opts_custom.fasta_ref)

        # Primer removal
        if (self._opts_custom.fasta_primers is not None) and (self._opts_custom.species != 'sars_cov_2'):
            raise ValueError("Primer removal is currently only supported for SARS-CoV-2")

    @staticmethod
    def _validate_ref_genome_file(path_fasta: Path) -> None:
        """
        Checks if the input file is a valid reference genome file.
        :param path_fasta: Input FASTA path
        :return: None
        """
        with path_fasta.open() as handle:
            for seq in SeqIO.parse(handle, 'fasta'):
                m = re.match(r'^[\w.]+-\w+', seq.id)
                if m:
                    continue
                raise ValueError(
                    f"Invalid reference genome, sequence '{seq.id}' does not match format ("
                    f"{{identifier}}-{{segment_name}}, use 'genome' as segment name for viral species without "
                    f"segments)")

    def determine_genome_size(self) -> int:
        """
        Determines the genome size.
        :return: Genome size
        """
        if self._opts_custom.fasta_ref is not None:
            logger.info('Extracting genome size from reference FASTA file')
            return fastautils.count_bases(self._opts_custom.fasta_ref)
        elif self._opts_custom.ref_genome_db is not None:
            logger.info('Extracting genome size from database')
            with open(self._opts_custom.ref_genome_db / 'genome_info.json') as handle:
                data_genome = json.load(handle)
                return data_genome['genome_size']
        raise ValueError('Cannot determine genome size')

    def _validate_config_data(self, config_data: dict) -> bool:
        """
        Validates the config data.
        :param config_data: Config data
        :return: True if valid, False otherwise
        """
        self.check_dbs(config_data)
        return True

    def _execute(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        self._validate_opts()

        # Build and validate the config file
        config_data = self._build_config()
        self._validate_config_data(config_data)

        # Create the config file and run snakefile
        self._script_out.dir.mkdir(parents=True, exist_ok=True)
        path_config = snakepipelineutils.generate_config_file(
            config_data, self._script_opts.working_dir
        )
        self.run_snakefile(path_config)

        # Copy the FASTA file of the consensus sequence (if specified)
        if self._script_out.fasta is not None:
            output_io_list = snakemakeutils.load_object(Path(
                self._script_opts.working_dir, iterativemapping.OUTPUT_FASTA_CONSENSUS_FINAL_TRIMMED))
            shutil.copyfile(output_io_list[0].path, self._script_out.fasta)

    def __config_add_yaml_data(self, config_data: dict) -> None:
        """
        Adds the data from parsing the config YAML file.
        :param config_data: Configuration data
        :return: None
        """
        path_config_template = Path(str(files('camel').joinpath('scripts/viralconsensuspipeline/config_data.yml')))
        logger.info(f'Adding config data from: {path_config_template}')
        with path_config_template.open() as handle_in:
            basepipeutils.dict_merge(
                config_data, yaml.safe_load(handle_in.read().format(
                    COV_MAX=self._script_opts.cov_max,
                    COV_MAX_SEGMENT=self._opts_custom.cov_max_segment,
                    DB_ROOT=config.dir_db,
                    GENOME_SIZE=self.determine_genome_size(),
                )))

    def __config_add_iterative_mapping_data(self, config_data: dict) -> None:
        """
        Adds the config data for the iterative mapping assay.
        :param config_data: Configuration data
        :return: None
        """
        # Clair3 model
        if self._opts_custom.clair3_model is not None:
            clair3_model = self._opts_custom.clair3_model
        elif self._script_in.type_ == model.InputType.ILLUMINA:
            logger.info("Clair3 model not specified, using default model for Illumina data")
            clair3_model = Path(config.dir_db, 'clair3', 'models', 'ilmn')
        else:
            logger.info("Clair3 model not specified, using default model for ONT data")
            clair3_model = Path(config.dir_db, 'clair3', 'models', 'ont')

        # Other values
        config_data['iterative_mapping'] = {
            'max_iter': self._opts_custom.max_iter,
            'variant_filters': {
                'min_af': self._opts_custom.variant_min_af,
                'min_dp': self._opts_custom.variant_min_dp,
                'min_qual': self._opts_custom.variant_min_qual,
                'min_mq': self._opts_custom.variant_min_mq},
            'clair3': {'model': str(clair3_model)}
        }

    def __config_add_coverage_data(self, config_data: dict) -> None:
        """
        Adds the config data for the downsampling & gap identification.
        :param config_data: Configuration data
        :return: None
        """
        config_data['low_depth'] = {
            'gap_depth_cutoff': self._opts_custom.gap_depth_cutoff,
            'gap_len_cutoff': self._opts_custom.gap_len_cutoff
        }

    def _build_config(self) -> dict:
        """
        Builds the configuration data for Snakemake.
        :return: Configuration data
        """
        # Input files and read type
        config_data = self.get_config_data() # TODO: check if invalid species raises an error

        # YAML data
        self.__config_add_yaml_data(config_data)
        config_data['analyses_selected'] = ['kraken2']
        config_data['analyses_selected'].extend(self._opts_custom.analyses)

        if self._opts_custom.species == 'other':
            config_data['contamination_check']['expected_species'] = self._opts_custom.species_name

        # Antiviral mutation detection (only for selected species)
        # species_antivirals = DATA_BY_SPECIES.get(self._opts_custom.species, {}).get('antivirals_species')
        # if species_antivirals is not None:
        #     config_data['analyses_selected'].append('antivirals')
        #     config_data['antivirals']['species'] = species_antivirals

        # Amplicon primer clipping
        if self._opts_custom.fasta_primers is not None:
            config_data['preprocess'] = {'ampligone': {'fasta': str(self._opts_custom.fasta_primers)}}
            config_data['analyses_selected'].append('ampligone')

        # Reference genome / reference selection
        config_data['fasta_ref'] = str(self._opts_custom.fasta_ref) if self._opts_custom.fasta_ref is not None else None
        config_data['ref_selection'] = {
            'db': str(self._opts_custom.ref_genome_db) if self._opts_custom.fasta_ref is None else None}

        # Other assays
        self.__config_add_coverage_data(config_data)
        self.__config_add_iterative_mapping_data(config_data)

        # Resolve species specific values
        config_data = basepipeutils.resolve_config(
            config_data, self._opts_custom.species
        )
        return config_data


@click.command(name='viral_consensus_pipeline', short_help='Extracts the consensus sequence from viral sequencing data')
@basescriptutils.add_input_opts()
@basescriptutils.add_output_opts
@basescriptutils.add_general_opts
@click.option('--analyses', type=str, help="Comma-separated list of analyses to run")
@cliutils.add_click_options_from_dataclass(Options, skip=['analyses'])
def main(**kwargs) -> None:
    """
    Extracts the consensus sequence from viral sequencing data.
    """
    script_input = basescriptutils.parse_script_input(kwargs)
    script_out = basescriptutils.parse_script_output(kwargs)
    script_opts = basescriptutils.parse_script_opts(kwargs)
    custom_opts = Options(
        analyses=kwargs['analyses'].split(',') if kwargs['analyses'] else [],
        **cliutils.from_kwargs(Options, kwargs, skip=['analyses'])
    )
    pipe_script = MainViralConsensusPipeline(script_input, script_out, script_opts, custom_opts)
    pipe_script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
