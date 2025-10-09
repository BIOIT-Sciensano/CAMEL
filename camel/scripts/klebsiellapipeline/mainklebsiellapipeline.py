#!/usr/bin/env python
import argparse
from collections.abc import Sequence
from typing import Optional

import yaml

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.klebsiellapipeline import CONFIG_DATA, SNAKEFILE_MAIN


class MainKlebsiellaPipeline(ReportPipeline):
    """
    Main class for running the Klebsiella pipeline.
    """

    CUSTOM_ANALYSES = [
        'kraken2', 'confindr', 'amrfinder', 'resfinder4', 'vfdb_core', 'kleborate', 'plasmidfinder', 'mob_suite',
        'bacmet', 'cgmlst', 'mlst', 'scgmlst', 'human_read_scrubbing', 'variant_calling']

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Klebsiella pipeline', '1.1', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return '<i>Klebsiella</i> pipeline'

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

    def __construct_config_file(self, input_files: dict[str, list[dict[str, str]]]) -> str:
        """
        Constructs the configuration file.
        :param input_files: Dictionary with the input files (keys can be FASTQ_PE, FASTQ_SE).
        :return: Configuration file
        """
        config_data = self.get_template_data(input_files)
        config_data['analyses'] = [key for key in MainKlebsiellaPipeline.CUSTOM_ANALYSES if vars(self._args)[key]]
        with CONFIG_DATA.open() as handle_in:
            mainscriptutils.dict_merge(config_data, yaml.load(handle_in.read().format(
                coverage_max=self._args.cov_max,
                export_bam='true' if self._args.report_include_bam else 'false',
                qc_typing_scheme='cgmlst' if self._args.cgmlst else 'mlst'
            ), Loader=yaml.SafeLoader))
        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        for analysis_key in MainKlebsiellaPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_', '-')}", action='store_true')
        return parser.parse_args(args)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainKlebsiellaPipeline()
    main.run()
