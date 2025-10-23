#!/usr/bin/env python
import argparse
from collections.abc import Sequence
from typing import Optional

import yaml

from camel.app.scriptutils import mainscriptutils
from camel.app.scriptutils.reportpipeline import ReportPipeline
from camel.app.core.snakemake import snakepipelineutils
from camel.app.loggers import initialize_logging
from camel.scripts.staphylococcuspipeline import CONFIG_DATA, SNAKEFILE_MAIN


class MainStaphylococcusPipeline(ReportPipeline):
    """
    Main class to run the Staphylococcus pipeline.
    """

    CUSTOM_ANALYSES = [
        'kraken2', 'confindr', 'rmlst', 'lrefinder', 'amrfinder', 'resfinder4', 'vfdb_core', 'virulencefinder', 'mlst',
        'cgmlst', 'spa_typing', 'sccmec_typing', 'plasmidfinder', 'mob_suite', 'se_toxins', 'bacmet',
        'human_read_scrubbing', 'variant_calling']

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Staphylococcus pipeline', '1.2', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return '<i>Staphylococcus</i> pipeline'

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
        config_data['analyses'] = [key for key in MainStaphylococcusPipeline.CUSTOM_ANALYSES if vars(self._args)[key]]

        with CONFIG_DATA.open() as handle_in:
            mainscriptutils.dict_merge(
                config_data, yaml.safe_load(handle_in.read().format(
                    coverage_max=self._args.cov_max,
                    qc_typing_scheme='cgmlst' if self._args.cgmlst else 'mlst',
                    export_bam='true' if self._args.report_include_bam else 'false',
                )))
        return snakepipelineutils.generate_config_file(config_data, self._args.working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        for analysis_key in MainStaphylococcusPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_', '-')}", action='store_true')
        return parser.parse_args(args)


if __name__ == '__main__':
    initialize_logging()
    main = MainStaphylococcusPipeline()
    main.run()
