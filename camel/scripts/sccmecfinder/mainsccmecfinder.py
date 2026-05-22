#!/usr/bin/env python
import dataclasses
import json
from pathlib import Path
from typing import Optional

import click
import yaml
from camelcore.app.reports.htmlreport import HtmlReport
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.utils import reportutils

from camel.app.cli import cliutils
from camel.app.loggers import initialize_logging, logger
from camel.app.scriptutils import inputhelper, model
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.app.scriptutils.inputhelper import helper_by_input_type
from camel.app.scriptutils.inputhelper.inputhelperbase import InputHelperBase
from camel.app.scriptutils.model import BaseOptions
from camel.app.wrappers.genedetectionwrapper import GeneDetectionWrapper


@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Custom options for the SCCmecFinder script.
    """
    db_mec_genes: Path = dataclasses.field(metadata={'help': 'Database containing mec genes'})
    profiles_mec_genes: Path = dataclasses.field(metadata={'help': 'Profiles for the mec genes'})
    working_dir: Path = dataclasses.field(default=Path.cwd(), metadata={'help': 'Working directory'})
    threads: int = dataclasses.field(default=4, metadata={'help': 'Number of threads to use'})


class MainSCCmecFinder(BaseScript[ScriptInput, ScriptOutput, Options]):
    """
    This tool is used to run the SCCmecFinder tool.
    """

    def __init__(self, script_in: ScriptInput, helper: InputHelperBase, script_out: ScriptOutput,
                 script_opts: Options) -> None:
        """
        Initializes the main script.
        :param script_in: Script input
        :param helper: Input helper
        :param script_out: Script output
        :param script_opts: Script options
        """
        super().__init__(
            name='SCCmecFinder',
            version='1.0.0',
            script_in=script_in,
            script_out=script_out,
            script_opts=script_opts
        )
        self._helper: InputHelperBase = helper
        self._report: HtmlReport | None = None
        self._informs: dict[str, str] = {}

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
        return None

    def _execute(self) -> None:
        """
        Executes the script.
        :return: None
        """
        # Init report
        self._report = reportutils.init_report(
            path_out=self._script_out.html,
            dir_out=self._script_out.dir,
            key='SCCmecFinder ouptput',
            title='SCCmecFinder (local)')
        self._report.add_html_object(reportutils.create_overview_section(
            version=self._version,
            dataset_name=self._script_in.name,
            input_file_str=self._script_in.input_str,
        ))
        self._report.save()

        # Run tools
        fasta_file = self._helper.prepare_fasta_input(self._script_in, self._report)
        detected_genes = self.__run_blast(fasta_file)

        # Save the output
        report_mec_type = self.__get_mec_type_overview(detected_genes)
        self._helper.export_output_and_commands_section(self._report, report_mec_type)

        # Save the JSON output (if specified)
        if self._script_out.json is not None:
            with self._script_out.json.open('w') as handle:
                json.dump(self._informs, handle, indent=2)
            logger.info(f'Informs exported to: {self._script_out.json}')

    def __run_blast(self, fasta_file: Path) -> list[str]:
        """
        Runs BLAST on the mec genes database.
        :param fasta_file: Input FASTA file
        :return: List of detected genes
        """
        wrapper = GeneDetectionWrapper(self._script_opts.working_dir / 'meca')
        wrapper.run_blast(
          fasta_file, self._script_in.name, {'path': str(self._script_opts.db_mec_genes)},
            self._script_opts.threads)
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
        with open(self._script_opts.profiles_mec_genes) as handle:
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

@click.command(name='sccmec_finder', short_help='Detection of SCCmec elements in Staphylococcus aureus')
@basescriptutils.add_input_opts(supported=[model.InputType.FASTA, model.InputType.ILLUMINA, model.InputType.ONT])
@basescriptutils.add_output_opts
@cliutils.add_click_options_from_dataclass(Options)
@inputhelper.add_helper_opts
def main(**kwargs) -> None:
    """
    Detection of SCCmec elements in Staphylococcus aureus.
    """
    # Parse the script input
    script_input = basescriptutils.parse_script_input(kwargs)
    script_opts = Options(**cliutils.from_kwargs(Options, kwargs))

    # Initialize the helper class to prepare the input
    helper = helper_by_input_type[script_input.type_](dir_=script_opts.working_dir, name=script_input.name)
    helper.set_opts(*helper.opts_from_cli(kwargs))

    # Run the main script
    script = MainSCCmecFinder(
        script_in=basescriptutils.parse_script_input(kwargs),
        script_out=basescriptutils.parse_script_output(kwargs),
        script_opts=script_opts,
        helper=helper
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
