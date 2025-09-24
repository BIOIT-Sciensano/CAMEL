#!/usr/bin/env python
import argparse
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

import pandas as pd
import pkg_resources
import yaml

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.galaxy.galaxyutils import GalaxyUtils
from camel.app.components.pipelines import absolute_path_by_pathlib
from camel.app.components.pipelines.basepipeline import BasePipeline
from camel.app.loggers import logger
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.hybridassemblypipeline import CONFIG_DATA


class MainHybridAssemblyPipeline(BasePipeline):
    """
    Main class to run the hybrid assembly pipeline.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the hybrid assembly pipeline class.
        :param args: arguments to be parsed
        :return: None
        """
        self._data_models = pd.read_table(pkg_resources.resource_filename(
            'camel', 'scripts/hybridassemblypipeline/basecalling_models.tsv'))
        self._args = MainHybridAssemblyPipeline._parse_arguments(args)
        _path_snakefile = pkg_resources.resource_filename('camel', 'scripts/hybridassemblypipeline/snakefile/main.smk')
        super().__init__(
            'Hybrid assembly pipeline', '0.1', _path_snakefile, args)
        if self._args.output_dir is None:
            self._args.output_dir = self._args.output_html.parent
        self._sample_name = MainHybridAssemblyPipeline._determine_sample_name(self._args)

    @staticmethod
    def _determine_sample_name(args: argparse.Namespace) -> str:
        """
        Determines the sample name from the provided arguments.
        :param args: Command line arguments
        :return: Parse sample name
        """
        if args.sample_name is not None:
            return args.sample_name
        logger.debug('Sample name not provided, determining name from short-reads')
        if args.fastq_pe_names is not None:
            return GalaxyUtils.determine_sample_name_from_fq([Path(fq) for fq in args.fastq_pe_names])
        return GalaxyUtils.determine_sample_name_from_fq([Path(fq) for fq in args.fastq_pe])

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        BasePipeline.add_common_arguments(argument_parser)

        # Output
        argument_parser.add_argument('--output-html', type=absolute_path_by_pathlib, required=True)
        argument_parser.add_argument('--output-dir', type=absolute_path_by_pathlib)

        # Parameters
        argument_parser.add_argument('--ploidy', type=int, choices=[1, 2], default=1)
        argument_parser.add_argument('--ont-qual', type=str, required=True,
                                     choices=['nano-corr', 'nano-hq', 'nano-raw'], default='nano-corr')
        argument_parser.add_argument('--flye-meta', action='store_true', default=None,
                                     help='enables the --meta option for the Flye assembly')
        argument_parser.add_argument('--expected-genome-size', type=str, default='')
        argument_parser.add_argument('--ont-basecalling-model', type=str, choices=[
            'r1041_e82_260bps_fast_g632', 'r1041_e82_260bps_hac_g632', 'r1041_e82_260bps_sup_g632',
            'r1041_e82_400bps_fast_g615', 'r1041_e82_400bps_fast_g632', 'r1041_e82_400bps_hac_g615',
            'r1041_e82_400bps_hac_g632', 'r1041_e82_400bps_sup_g615', 'r104_e81_hac_g5015', 'r104_e81_sup_g5015',
            'r104_e81_sup_g610', 'r1041_e82_400bps_sup_v5.0.0', 'r941_min_hac_g507', 'r941_min_high_g360',
            'r941_min_sup_g507', 'r941_prom_hac_g507', 'r941_prom_high_g360', 'r941_prom_sup_g507'],
                                     default='r941_prom_sup_g507')
        argument_parser.add_argument('--filtlong-keep-percent', type=int)

        # Variant calling
        argument_parser.add_argument('--freebayes-min-alternate-fraction', type=float, default=0.5)
        argument_parser.add_argument('--freebayes-min-alternate-count', type=int, default=10)
        argument_parser.add_argument('--clair3-long-indel', action='store_true', default=None)
        argument_parser.add_argument('--sniffles-min-svlen', type=int)
        argument_parser.add_argument('--sniffles-mapq', type=int)
        argument_parser.add_argument('--sniffles-min-support', type=int)

        # Flag to enable Unicycler
        argument_parser.add_argument('--unicycler', action='store_true', default=False)

        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the hybrid assembly pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        path_config = self.__create_snakemake_config_data(input_files)
        self._run_snakemake_main(path_config)

    def __create_snakemake_config_data(self, input_files: dict[str, list[dict[str, str]]]) -> str:
        """
        Creates a Snakemake configuration file.
        :input_files: Dictionary with the input files (keys can be FASTQ_PE, FASTQ_SE).
        :return: Config file data
        """
        config_data = self.get_template_data(input_files)

        # Add report-specific entries
        mainscriptutils.dict_merge(config_data, {
            'output_dir': str(self._args.output_dir.absolute()),
            'output_html': str(self._args.output_html.absolute())
        })

        # Assembly steps
        config_data['assembly_steps'] = ['Flye', 'Medaka', 'Polypolish', 'POLCA']
        config_data['base_assemblies'] = ['flye']
        if self._args.unicycler:
            config_data['assembly_steps'].extend(['Unicycler', 'Medaka-Unicycler', 'Polypolish-Unicycler', 'POLCA-Unicycler'])
            config_data['base_assemblies'].extend(['unicycler'])

        with open(CONFIG_DATA) as handle_in:
            mainscriptutils.dict_merge(
                config_data, yaml.safe_load(handle_in.read().format(
                    clair3_haploid_precise=True if self._args.ploidy == 1 else False,
                    clair3_long_indel=True if self._args.clair3_long_indel is not None else False,
                    clair3_model_path=str(self._get_clair3_model(self._args.ont_basecalling_model)),
                    freebayes_min_alternate_count=self._args.freebayes_min_alternate_count,
                    freebayes_min_alternate_fraction=self._args.freebayes_min_alternate_fraction,
                    freebayes_ploidy=self._args.ploidy,
                    expected_size=self._args.expected_genome_size if self._args.expected_genome_size else False,
                    is_meta=True if self._args.flye_meta else False,
                    medaka_model=self._args.ont_basecalling_model,
                    nano_corr=True if self._args.ont_qual == 'nano-corr' else False,
                    nano_hq=True if self._args.ont_qual == 'nano-hq' else False,
                    nano_raw=True if self._args.ont_qual == 'nano-raw' else False,
                    sniffles_mapq=self._args.sniffles_mapq if self._args.sniffles_mapq is not None else 40,
                    sniffles_min_support=self._args.sniffles_min_support if self._args.sniffles_min_support is not None else 'auto',
                    sniffles_min_svlen=self._args.sniffles_min_svlen if self._args.sniffles_min_svlen is not None else 35
                )))
        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    def _get_clair3_model(self, medaka_model: str) -> Path:
        """
        Gets the most suitable clair3 model for the given medaka model.
        From documentation: the medaka model with the highest version equal to or less than the guppy version should
        be selected.
        :param medaka_model: Medaka basecalling model
        :return: Path to the Clair3 model
        """
        clair3_model = next(
            r['clair3_model'] for r in self._data_models.to_dict('records') if r['medaka_model'] == medaka_model)
        logger.info(f"Best matching Clair3 model for medaka model '{medaka_model}': {clair3_model}")
        path_model = Path(Camel.get_instance().config['db_root'], 'clair3', 'models', clair3_model)
        if not path_model.exists():
            raise ValueError(f'Clair3 model not found: {path_model}')
        return path_model


if __name__ == '__main__':
    Camel.get_instance()
    main = MainHybridAssemblyPipeline()
    main.run()
