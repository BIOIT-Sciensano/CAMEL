#!/usr/bin/env python
import argparse
from collections.abc import Sequence
from typing import Optional

import yaml

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.salmonellapipeline import CONFIG_DATA, SNAKEFILE_MAIN


class MainSalmonellaPipeline(ReportPipeline):
    """
    Main class to run the Salmonella pipeline.
    """

    CUSTOM_ANALYSES = [
        'kraken2', 'confindr', 'mlst', 'cgmlst', 'rmlst', 'resfinder4', 'serotype', 'vfdb_core', 'spifinder',
        'mykrobe', 'abritamr', 'mob_suite', 'human_read_scrubbing', 'variant_calling', 'amrfinder']

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Salmonella pipeline', '1.0', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return '<i>Salmonella</i> pipeline'

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
        config_data['analyses'] = [key for key in MainSalmonellaPipeline.CUSTOM_ANALYSES if vars(self._args)[key]]
        with CONFIG_DATA.open() as handle_in:
            mainscriptutils.dict_merge(
                config_data, yaml.safe_load(handle_in.read().format(
                qc_typing_scheme='cgmlst' if self._args.cgmlst else 'mlst',
                export_bam='true' if self._args.report_include_bam else 'false',
                coverage_max=self._args.cov_max
            )))
        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        for analysis_key in MainSalmonellaPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_', '-')}", action='store_true')
        return parser.parse_args(args)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainSalmonellaPipeline()
    main.run()
