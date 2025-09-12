#!/usr/bin/env python
import argparse
import json
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

import yaml

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.workflows.genedetectionwrapper import GeneDetectionWrapper
from camel.app.components.workflows.inputtype import helper_by_input_type
from camel.app.loggers import logger


class MainSCCmecFinder:
    """
    This tool is used to run the SCCmecFinder tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        :param args: Arguments, if not set they are removed from the command line
        """
        self._camel = Camel()
        self._args = MainSCCmecFinder._parse_arguments(args)
        self._sample_name = mainscriptutils.determine_sample_name(self._args)
        self._helper = helper_by_input_type[self._args.input_type](self._args.working_dir, self._sample_name)
        self._report = None
        self._informs = {}

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        mainscriptutils.add_common_arguments(argument_parser)
        mainscriptutils.add_input_files_arguments(argument_parser)
        mainscriptutils.add_assembly_arguments(argument_parser)
        argument_parser.add_argument('--db-mec-genes', type=Path, help="Database containing mec genes.", required=True)
        argument_parser.add_argument('--profiles-mec-genes', help="Profiles for the mec genes", required=True)
        argument_parser.add_argument('--output-json', type=Path, help='Output path to store informs')
        return argument_parser.parse_args(args)

    @staticmethod
    def __get_matching_complex(detected_genes: list[str], genes_by_complex: dict[str, list[str]]) -> \
            Optional[str]:
        """
        Returns the matching complex (if there is one).
        :param genes_by_complex: Genes by complex
        :return: Complex (or None if there is none found)
        """
        for complex_, genes in genes_by_complex.items():
            if all(g in detected_genes for g in genes):
                return complex_
        logger.debug("No complex found")

    def run(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Init report
        self._report = mainscriptutils.init_report(
            self._args.output_html, self._args.output_dir, 'SCCmecFinder ouptput', 'SCCmecFinder (local)')
        self._report.add_html_object(mainscriptutils.generate_analysis_info_section(self._args))
        self._report.save()

        # Run tools
        fasta_file = self._helper.prepare_fasta_input(self._report, self._args)
        detected_genes = self.__run_blast(fasta_file)

        # Save the output
        report_mec_type = self.__get_mec_type_overview(detected_genes)
        self._helper.export_output_and_commands_section(self._report, report_mec_type)

        # Save the JSON output (if specified)
        if self._args.output_json is not None:
            with self._args.output_json.open('w') as handle:
                json.dump(self._informs, handle, indent=2)
            logger.info(f'Informs exported to: {self._args.output_json}')

    def __run_blast(self, fasta_file: Path) -> list[str]:
        """
        Runs BLAST on the mec genes database.
        :param fasta_file: Input FASTA file
        :return: List of detected genes
        """
        wrapper = GeneDetectionWrapper(self._args.working_dir / 'meca')
        wrapper.run_workflow_blast(
          fasta_file, self._sample_name, {'path': str(self._args.db_mec_genes)}, self._args.threads)
        self._report.add_html_object(wrapper.output.report_section)
        wrapper.output.report_section.copy_files(self._report.output_dir)
        self._helper.informs.append(wrapper.output.informs)
        self._report.save()
        return [d.locus.split(':')[0] for d in wrapper.output.detected_hits]

    def __get_mec_type_overview(self, detected_genes: list[str]) -> HtmlReportSection:
        """
        Determines the mec type based on the detected genes and adds it to the report.
        :param detected_genes: Detected genes
        :return: None
        """
        with open(self._args.profiles_mec_genes) as handle:
            profiles = yaml.safe_load(handle)
        self._informs['ccr_complex'] = MainSCCmecFinder.__get_matching_complex(
            detected_genes, profiles['ccr_genes_complexes'])
        self._informs['mec_complex'] = MainSCCmecFinder.__get_matching_complex(
            detected_genes, profiles['mec_genes_complexes'])
        self._informs['sccmec_type'] = MainSCCmecFinder.__get_matching_complex(
            detected_genes, profiles['SCC_mec_types'])
        section = HtmlReportSection('SCCmec type', 3)
        section.add_table([
            ['SCCmec type:', self._informs['sccmec_type'] if self._informs['sccmec_type'] is not None else '-'],
            ['<i>mec</i> class:', self._informs['mec_complex'] if self._informs['mec_complex'] is not None else '-'],
            ['<i>ccr</i> class:', self._informs['sccmec_type'] if self._informs['sccmec_type'] is not None else '-']
        ], None, [('class', 'information')])
        return section


if __name__ == '__main__':
    Camel.get_instance()
    main = MainSCCmecFinder()
    main.run()
