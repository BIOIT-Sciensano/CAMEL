#!/usr/bin/env python
import argparse
from collections.abc import Sequence
from typing import Optional

import yaml

from camel.app.scriptutils.reportpipeline import ReportPipeline
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.snakemake import snakepipelineutils
from camel.app.loggers import initialize_logging
from camel.scripts.abritamr import CONFIG_DATA, SNAKEFILE_MAIN


class MainAbriTAMR(ReportPipeline):
    """
    Main class to run the AbriTAMR standalone pipeline.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        :return: None
        """
        super().__init__('AbriTAMR standalone', '0.2', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return 'AbriTAMR'

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        if self._args.input_type != 'fasta':
            raise InvalidToolInputError('The input type for this tool has to be fasta')
        input_files = self._symlink_input()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)

    def __construct_config_file(self, input_files: dict[str, list[dict[str, str]]]) -> str:
        """
        Constructs the configuration file.
        :param input_files: Dictionary with the input files (key can only be FASTA).
        :return: Configuration file
        """
        config_data = self.get_template_data(input_files)
        # Add existing config data
        with open(CONFIG_DATA) as handle_in:
            config_data.update(yaml.safe_load(handle_in.read()))
        config_data['abritamr']['species'] = self._args.species
        config_data['analyses'] = ['abritamr']
        return snakepipelineutils.generate_config_file(config_data, self._args.working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments (optional)
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(argument_parser)
        argument_parser.add_argument('--species', required=True, choices=[
            'Acinetobacter_baumannii', 'Burkholderia_cepacia', 'Burkholderia_pseudomallei', 'Burkholderia_mallei',
            'Campylobacter', 'Citrobacter_freundii', 'Clostridioides_difficile', 'Corynebacterium_diphtheriae',
            'Enterobacter_asburiae', 'Enterobacter_cloacae', 'Enterococcus_faecalis', 'Enterococcus_faecium',
            'Escherichia', 'Klebsiella_oxytoca', 'Klebsiella_pneumoniae', 'Neisseria_gonorrhoeae',
            'Neisseria_meningitidis', 'Pseudomonas_aeruginosa', 'Salmonella', 'Serratia_marcescens',
            'Staphylococcus_aureus', 'Staphylococcus_pseudintermedius', 'Streptococcus_agalactiae',
            'Streptococcus_pneumoniae', 'Streptococcus_pyogenes', 'Vibrio_cholerae', 'Vibrio_vulfinicus',
            'Vibrio_parahaemolyticus'])
        return argument_parser.parse_args(args)


if __name__ == '__main__':
    initialize_logging()
    main = MainAbriTAMR()
    main.run()
