#!/usr/bin/env python
import argparse
from collections.abc import Sequence
from typing import Optional

import yaml

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.loggers import logger
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.bacilluspipeline import CONFIG_DATA, SNAKEFILE_MAIN


class MainBacillusPipeline(ReportPipeline):
    """
    Main class to run the Bacillus pipeline.
    """

    CUSTOM_ANALYSES = {
        'common': ['rmlst', 'plasmidfinder', 'mobsuite', 'vfdb_core', 'amrfinder', 'kraken2', 'confindr',
                   'human_read_scrubbing', 'variant_calling'],
        'cereus': ['btyper', 'mlst_cereus', 'cgmlst_cereus'],
        'subtilis': ['fastani', 'mlst_subtilis', 'gmo', 'straingst']
    }

    DATA_BY_SPECIES = {
        'cereus': {
            'full_name': 'Bacillus cereus',
            'ref_name': 'NZ_CP017060.1',
            'ref_fasta': '/db/refgenomes/Bacillus/NZ_CP017060.1.fasta',
            'ref_gff3': '/db/refgenomes/Bacillus/NZ_CP017060.1.gff3',
            'ref_url': 'https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP017060.1'
        },
        'subtilis': {
            'full_name': 'Bacillus subtilis',
            'ref_name': 'NC_000964.3',
            'ref_fasta': '/db/refgenomes/Bacillus/NC_000964.3.fasta',
            'ref_gff3': '/db/refgenomes/Bacillus/NC_000964.3.gff3',
            'ref_url': 'https://www.ncbi.nlm.nih.gov/nuccore/NC_000964.3'
        }
    }

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Bacillus pipeline', '0.1', SNAKEFILE_MAIN, args)
        self._args = MainBacillusPipeline._parse_arguments(args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return '<i>Bacillus</i> pipeline'

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        self._validate_input_files()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)
        self._export_assembly()

    def __get_qc_typing_scheme(self) -> str:
        """
        Returns the typing scheme used for QC.
        :return: Typing scheme
        """
        if self._args.species == 'cereus':
            return 'cgmlst_cereus' if self._args.cgmlst_cereus is True else 'mlst_cereus'
        else:
            return 'mlst_subtilis'

    def __construct_config_file(self, input_files: dict[str, list[dict[str, str]]]) -> str:
        """
        Constructs the configuration file.
        :input_files: Dictionary with the input files (keys can be FASTQ_PE, FASTQ_SE).
        :return: Configuration file
        """
        config_data = self.get_template_data(input_files)

        # Analyses to perform
        config_data['analyses'] = []
        for group, keys in MainBacillusPipeline.CUSTOM_ANALYSES.items():
            for key in keys:
                if not vars(self._args)[key]:
                    continue
                if group not in ('common', self._args.species):
                    logger.warning(f"Analysis '{key}' not supported for species '{self._args.species}'")
                    continue
                if key == 'variant_calling' and self._args.input_type not in ('illumina', 'fasta'):
                    continue
                if key == 'confindr' and self._args.input_type not in ('illumina', 'ont', 'hybrid'):
                    continue
                config_data['analyses'].append(key)

        # Parse template
        with open(CONFIG_DATA) as handle_in:
            mainscriptutils.dict_merge(config_data, yaml.load(handle_in.read().format(
                coverage_max=self._args.cov_max,
                export_bam='true' if self._args.report_include_bam else 'false',
                export_fastq='true' if self._args.report_include_fastq else 'false',
                mobsuite_contig_report=self._args.mobsuite_contig_report,
                qc_typing_scheme=self.__get_qc_typing_scheme(),
                ref_fasta=MainBacillusPipeline.DATA_BY_SPECIES[self._args.species].get('ref_fasta', 'null'),
                ref_gff=MainBacillusPipeline.DATA_BY_SPECIES[self._args.species].get('ref_gff', 'null'),
                ref_name=MainBacillusPipeline.DATA_BY_SPECIES[self._args.species].get('ref_name', 'null'),
                ref_url=MainBacillusPipeline.DATA_BY_SPECIES[self._args.species].get('ref_url', 'null'),
                species=self._args.species,
            ), Loader=yaml.SafeLoader))

        # ONT settings
        if self._args.input_type in ['ont', 'hybrid']:
            config_data['assembly']['flye'] = {
                **config_data['assembly'].get('flye', {})}

        # Illumina settings
        if self._args.library is not None:
            config_data['read_trimming']['adapter'] = self._args.library

        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: command line arguments
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        parser.add_argument('--species', type=str, choices=['cereus', 'subtilis'], required=True)
        parser.add_argument('--mobsuite-contig-report', action='store_true')
        for _, keys in MainBacillusPipeline.CUSTOM_ANALYSES.items():
            for analysis_key in keys:
                parser.add_argument(f"--{analysis_key.replace('_', '-')}", action='store_true')
        return parser.parse_args(args)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainBacillusPipeline()
    main.run()
