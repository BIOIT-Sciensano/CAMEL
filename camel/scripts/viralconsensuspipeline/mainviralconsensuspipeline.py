#!/usr/bin/env python
import argparse
import json
import re
import shutil
from collections.abc import Sequence
from importlib.resources import files
from pathlib import Path
from typing import Optional

import yaml
from Bio import SeqIO

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.loggers import logger
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.viralconsensuspipeline import SNAKEFILE_MAIN
from camel.scripts.viralconsensuspipeline.snakefile import iterativemapping


class MainViralConsensusPipeline(ReportPipeline):
    """
    Main script for the viral consensus pipeline.
    """

    CUSTOM_ANALYSES = []
    DB_ROOT = Path(Camel.get_instance().config['db_root'], 'pipelines', 'viral_consensus', 'version_1.1')

    SUPPORTED_SPECIES = {
        'influenza_a': {
            'antivirals_species': 'A',
            'name': 'Influenza A',
            'k2_name': 'Alphainfluenzavirus influenzae',
            'nextclade_mash_db': str(DB_ROOT / 'subtype_mash' / 'influenza_a'),
            'nextclade_capitalize': True,
        },
        'influenza_b': {
            'antivirals_species': 'B',
            'name': 'Influenza B',
            'k2_name': 'Betainfluenzavirus influenzae',
            'nextclade_segments': [],
            'nextclade_mash_db': str(DB_ROOT / 'subtype_mash' / 'influenza_b'),
            'nextclade_capitalize': True,
        },
        'sars_cov_2': {
            'antivirals_species': None,
            'name': 'SARS-CoV-2',
            'k2_name': 'Severe acute respiratory syndrome-related coronavirus',
            'nextclade_dbs': {
                'genome': str(Path(Camel.get_instance().config['db_root'], 'nextclade3', 'sars-cov-2',))
            },
        },
        'rsv_a': {
            'name': 'Respiratory syncytial virus (A)',
            'k2_name': 'Orthopneumovirus hominis',
            'nextclade_dbs': {
                'genome': str(Path(Camel.get_instance().config['db_root'], 'nextclade3', 'rsv_a'))
            }
        },
        'rsv_b': {
            'name': 'Respiratory syncytial virus (B)',
            'k2_name': 'Orthopneumovirus hominis',
            'nextclade_dbs': {
                'genome': str(Path(Camel.get_instance().config['db_root'], 'nextclade3', 'rsv_b'))
            }
        },
    }

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        """
        super().__init__('Viral consensus pipeline', '1.1', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return 'Viral consensus pipeline'

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Command line arguments
        :return: Parsed arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)

        # Species
        parser.add_argument('--species', required=True)
        parser.add_argument('--species-name', help='Species name (used for Kraken 2 classification)')

        # Reference genome & db
        parser.add_argument('--fasta-ref', type=Path, help='(Optional) reference genome')
        parser.add_argument('--fasta-ref-name', help='Name of the reference genome FASTA file')
        parser.add_argument('--ref-genome-db', type=Path,
                            help='Database with reference genomes (for automatic reference genome selection)')
        # Human read scrubbing
        parser.add_argument('--human-read-scrubbing', action='store_true', help='Remove human reads at the start of the pipeline')
        # Primer removal
        parser.add_argument('--fasta-primers', type=Path, help='Path to FASTA file with primer sequences')
        parser.add_argument('--fasta-primers-name', type=str, help='Name of the FASTA file with primer sequences')
        # Downsampling & gaps
        parser.set_defaults(cov_max=100_000)
        parser.add_argument('--cov-max-segment', type=int, default=10_000, help='Maximum coverage (per segment)')
        parser.add_argument(
            '--gap-depth-cutoff', type=int, default=5,
            help='Positions with a depth smaller than this value are flagged as missing / gaps')
        parser.add_argument(
            '--gap-len-cutoff', type=int, default=10, help='Minimum length to mark a region as a gap')

        # Iterative mapping & variant filtering
        parser.add_argument('--max-iter', type=int, default=6, help='Maximum number of iterations')
        parser.add_argument('--variant-min-af', type=float, default=0.5, help='Minimum allele frequency')
        parser.add_argument('--variant-min-dp', type=float, default=5, help='Minimum depth at variant position')
        parser.add_argument('--variant-min-qual', type=int, default=10, help='Minimum variant quality')
        parser.add_argument('--variant-min-mq', type=int, default=30, help='Minimum mapping quality')
        parser.add_argument('--clair3-model', type=Path, help='Clair3 variant calling model')

        # Custom analyses
        for analysis_key in MainViralConsensusPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_', '-')}", action='store_true')

        return parser.parse_args(args)

    @staticmethod
    def _check_args(args: argparse.Namespace) -> None:
        """
        Checks if the provided arguments are valid.
        :param args: Command line-arguments
        :return: None
        """
        # Supported species
        supported_species = list(MainViralConsensusPipeline.SUPPORTED_SPECIES.keys())
        if (args.species not in supported_species) and (args.species != 'other'):
            raise ValueError(f"Unsupported species: {args.species} (options: {','.join(supported_species)})")

        # Species name
        if (args.species == 'other') and (args.species_name is None):
            raise ValueError("Species name needs to be specified when species is 'other'")

        # Reference genome
        if args.fasta_ref is not None:
            MainViralConsensusPipeline._validate_ref_genome_file(args.fasta_ref)

        # Primer removal
        if (args.fasta_primers is not None) and (args.species != 'sars_cov_2'):
            raise ValueError('Primer removal is currently only supported for SARS-CoV-2')

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
        if self._args.fasta_ref is not None:
            logger.info('Extracting genome size from reference FASTA file')
            with open(self._args.fasta_ref) as handle:
                return sum(len(s) for s in SeqIO.parse(handle, 'fasta'))
        elif self._args.ref_genome_db is not None:
            logger.info('Extracting genome size from database')
            with open(self._args.ref_genome_db / 'genome_info.json') as handle:
                data_genome = json.load(handle)
                return data_genome['genome_size']
        raise ValueError('Cannot determine genome size')

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        # Check arguments
        MainViralConsensusPipeline._check_args(self._args)

        # Create symlinks for the input files
        if self._args.fasta_ref_name is not None:
            path_fasta = self._args.working_dir / self._args.fasta_ref_name
            path_fasta.symlink_to(self._args.fasta_ref)
            self._args.fasta_ref = path_fasta
        if self._args.fasta_primers_name is not None:
            path_fasta_primers = self._args.working_dir / self._args.fasta_primers_name
            path_fasta_primers.symlink_to(self._args.fasta_primers)
            self._args.fasta_primers = path_fasta_primers
        input_files = self._symlink_input()

        # Create config file and run snakemake
        path_config = self.__construct_config_file(input_files)
        self._run_snakemake_main(str(path_config))

        # Copy the FASTA file of the consensus sequence (if specified)
        if self._args.output_fasta is not None:
            output_io_list = SnakemakeUtils.load_object(Path(
                self._args.working_dir, iterativemapping.OUTPUT_ITERATIVE_MAPPING_FASTA_CONSENSUS_FINAL_TRIMMED))
            shutil.copyfile(output_io_list[0].path, self._args.output_fasta)

    def __config_add_yaml_data(self, config_data: dict) -> None:
        """
        Adds the data from parsing the config YAML file.
        :param config_data: Configuration data
        :return: None
        """
        path_config_template = Path(str(files('camel').joinpath('scripts/viralconsensuspipeline/config_data.yml')))
        logger.info(f'Adding config data from: {path_config_template}')
        with path_config_template.open() as handle_in:
            mainscriptutils.dict_merge(
                config_data, yaml.safe_load(handle_in.read().format(
                    export_fastq='true' if self._args.report_include_fastq else 'false',
                    cov_max=self._args.cov_max,
                    cov_max_segment=self._args.cov_max_segment,
                    genome_size=self.determine_genome_size(),
                    expected_species=MainViralConsensusPipeline.SUPPORTED_SPECIES[self._args.species]['k2_name'] if
                    self._args.species != 'other' else self._args.species_name
                )))

    def __config_add_nextclade_data(self, config_data: dict) -> None:
        """
        Adds the config data for the nextclade assay.
        :param config_data: Configuration data
        :return: None
        """
        data_by_species = MainViralConsensusPipeline.SUPPORTED_SPECIES
        config_data['nextclade'] = {}
        if data_by_species.get(self._args.species, {}).get('nextclade_dbs') is not None:
            config_data['nextclade']['dbs'] = data_by_species[self._args.species]['nextclade_dbs']
            config_data['analyses'].append('nextclade')
        elif data_by_species.get(self._args.species, {}).get('nextclade_mash_db') is not None:
            config_data['nextclade']['db_mash'] = str(data_by_species[self._args.species]['nextclade_mash_db'])
            config_data['analyses'].append('nextclade')
        if data_by_species.get(self._args.species, {}).get('nextclade_capitalize', False) is True:
            config_data['nextclade']['capitalize'] = True

    def __config_add_iterative_mapping_data(self, config_data: dict) -> None:
        """
        Adds the config data for the iterative mapping assay.
        :param config_data: Configuration data
        :return: None
        """
        # Clair3 model
        if self._args.clair3_model is not None:
            clair3_model = self._args.clair3_model
        elif self._args.input_type == 'illumina':
            logger.info("Clair3 model not specified, using default model for Illumina data")
            clair3_model = Path(Camel.get_instance().config['db_root'], 'clair3', 'models', 'ilmn')
        else:
            logger.info("Clair3 model not specified, using default model for ONT data")
            clair3_model = Path(Camel.get_instance().config['db_root'], 'clair3', 'models', 'ont')

        # Other values
        config_data['iterative_mapping'] = {
            'max_iter': self._args.max_iter,
            'variant_filters': {
                'min_af': self._args.variant_min_af,
                'min_dp': self._args.variant_min_dp,
                'min_qual': self._args.variant_min_qual,
                'min_mq': self._args.variant_min_mq},
            'clair3': {'model': str(clair3_model)}
        }

    def __config_add_coverage_data(self, config_data: dict) -> None:
        """
        Adds the config data for the downsampling & gap identification.
        :param config_data: Configuration data
        :return: None
        """
        config_data['low_depth'] = {
            'gap_depth_cutoff': self._args.gap_depth_cutoff,
            'gap_len_cutoff': self._args.gap_len_cutoff
        }

    def __construct_config_file(self, input_files: dict[str, list[dict[str, str]]]) -> Path:
        """
        Constructs the Snakemake configuration file.
        :param input_files: Input FASTQ files
        :return: Path to the config file
        """
        # Input files & read type
        config_data = self.get_template_data(input_files)

        # YAML data
        self.__config_add_yaml_data(config_data)
        config_data['analyses'] = ['kraken2']
        if self._args.human_read_scrubbing:
            config_data['analyses'].append('human_read_scrubbing')

        # Antiviral mutation detection (only for selected species)
        species_antivirals = MainViralConsensusPipeline.SUPPORTED_SPECIES.get(self._args.species, {}).get('antivirals_species')
        if species_antivirals is not None:
            config_data['analyses'].append('antivirals')
            config_data['antivirals'] = {'species': species_antivirals}

        # Amplicon primer clipping
        if self._args.fasta_primers is not None:
            config_data['preprocess'] = {'ampligone': {'fasta': str(self._args.fasta_primers)}}
            config_data['analyses'].append('ampligone')

        # Reference genome / reference selection
        config_data['fasta_ref'] = str(self._args.fasta_ref) if self._args.fasta_ref is not None else None
        config_data['ref_selection'] = {'db': str(self._args.ref_genome_db) if self._args.fasta_ref is None else None}

        # Other assays
        self.__config_add_coverage_data(config_data)
        self.__config_add_nextclade_data(config_data)
        self.__config_add_iterative_mapping_data(config_data)

        # Create YAML file
        path_config = SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)
        return Path(path_config)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainViralConsensusPipeline()
    main.run()
