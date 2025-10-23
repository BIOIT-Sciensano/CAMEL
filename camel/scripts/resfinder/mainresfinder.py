#!/usr/bin/env python
import argparse
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

from camel.app.core.reports import reportutils
from camel.app.scriptutils import mainscriptutils
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import initialize_logging
from camel.app.tools.resfinder.resfinder import ResFinder
from camel.app.tools.resfinder.resfinderreporter import ResFinderReporter


class MainResFinder:
    """
    This class is used to run the main ResFinder local script.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        :param args: Arguments (optional)
        """
        self._args = MainResFinder.parse_arguments(args)
        self._sample_name = mainscriptutils.determine_sample_name(self._args)

    @staticmethod
    def parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()

        mainscriptutils.add_common_arguments(argument_parser)

        argument_parser.add_argument('--fasta', help="Input FASTA file", type=Path)
        argument_parser.add_argument('--fasta-name', help="Input FASTA file name", type=str)

        argument_parser.add_argument('--fastq-pe', help="Input PE FASTQ files", nargs=2, type=Path)
        argument_parser.add_argument('--fastq-pe-names', help="Input PE FASTQ file names", nargs=2, type=Path)

        argument_parser.add_argument('--fastq-se', help="Input SE FASTQ file", type=Path)
        argument_parser.add_argument('--fastq-se-name', help="Input SE FASTQ file name")

        argument_parser.add_argument('--db-directory', help="Path containing the resfinder and pointfinder dbs",
                                     type=Path, required=True)

        argument_parser.add_argument('--point', action='store_true', default=None)
        argument_parser.add_argument('--acquired', action='store_true', default=None)

        argument_parser.add_argument('--min-cov', type=int, default=60,
                                     help='Minimum (breadth-of) coverage of ResFinder')
        argument_parser.add_argument('--threshold', type=int, default=80,
                                     help='Threshold for identity of ResFinder')
        argument_parser.add_argument('--acq-overlap', type=int, default=30,
                                     help=' Genes are allowed to overlap this number of nucleotides. Default: 30.')

        argument_parser.add_argument('--species', choices=[
            'Campylobacter', 'Campylobacter_jejuni', 'Campylobacter_coli', 'Enterococcus_faecalis',
            'Enterococcus_faecium', 'Escherichia_coli', 'Helicobacter_pylori', 'Klebsiella',
            'Mycobacterium_tuberculosis', 'Neisseria_gonorrhoeae', 'Plasmodium_falciparum', 'Salmonella',
            'Salmonella_enterica', 'Staphylococcus_aureus'], required=False, default=None)

        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        # Initialize report
        report = mainscriptutils.init_report(
            self._args.output_html, self._args.output_dir, 'ResFinder report', 'ResFinder')
        additional_info = [
            ['Species:', '<i>{}</i>'.format(self._args.species.replace('"', '')) if
                self._args.species is not None else 'Not specified'],
            ['Min % identity:',
             f'{self._args.threshold}' if self._args.threshold is not None else 'Curated (default)'],
            ['Min % coverage:', f'{self._args.min_cov}'],
        ]
        report.add_html_object(mainscriptutils.generate_analysis_info_section(self._args, additional_info))
        report.save()

        # Run tools
        resfinder = self.__run_resfinder()
        section = self.__run_reporter(resfinder)
        report.add_html_object(section)
        section.copy_files(report.output_dir)

        # Save report
        all_informs = [resfinder.informs]
        report.add_html_object(reportutils.create_commands_section(all_informs, self._args.working_dir))
        report.add_html_object(reportutils.create_citations_section(['Bortolaia_2020-resfinder_4.0']))
        report.save()

    def __run_resfinder(self) -> ResFinder:
        """
        Runs ResFinder.
        :return: ResFinder tool instance.
        """
        resfinder = ResFinder()
        resfinder.add_input_files({'DIR': [ToolIODirectory(self._args.db_directory)]})
        if self._args.fasta is not None:
            resfinder.add_input_files({'FASTA': [ToolIOFile(self._args.fasta)]})
        elif self._args.fastq_pe is not None:
            resfinder.add_input_files({'FASTQ_PE': [ToolIOFile(self._args.fastq_pe[0]),
                                                    ToolIOFile(self._args.fastq_pe[1])]})
        elif self._args.fastq_se is not None:
            resfinder.add_input_files({'FASTQ_SE': [ToolIOFile(self._args.fastq_se[0])]})

        resfinder.update_parameters(output_path=self._args.working_dir, min_cov=0.6, threshold=0.8)

        if self._args.min_cov != 60:
            resfinder.update_parameters(min_cov=self._args.min_cov / 100.0)
        if self._args.threshold != 80:
            resfinder.update_parameters(threshold=self._args.threshold / 100.0)
        if self._args.point is not None:
            try:
                resfinder.update_parameters(point=True, species='"' + self._args.species.replace('_', ' ') + '"')
            except AttributeError:
                raise ValueError('--point requires a --species argument')
        if self._args.acquired is not None:
            resfinder.update_parameters(acquired=True, acq_overlap=self._args.acq_overlap)

        resfinder.run(self._args.working_dir)
        return resfinder

    def __run_reporter(self, resfinder: ResFinder) -> HtmlReportSection:
        """
        Runs resfinder reporter.
        :param resfinder: ResFinder tool instance.
        :return: None.
        """
        reporter = ResFinderReporter()
        reporter.add_input_files({'TSV_pheno_general': resfinder.tool_outputs['TSV_pheno_general']})
        if self._args.acquired is not None:
            reporter.add_input_files({'TSV_genes': resfinder.tool_outputs['TSV_genes']})
        if self._args.point is not None:
            reporter.add_input_files({'TSV_point': resfinder.tool_outputs['TSV_point'],
                                      'TSV_pheno_species': resfinder.tool_outputs['TSV_pheno_species']})
        reporter.add_input_informs({'resfinder': resfinder.informs})
        reporter.run()
        return reporter.tool_outputs['VAL_HTML'][0].value


if __name__ == '__main__':
    initialize_logging()
    resfinder_main = MainResFinder()
    resfinder_main.run()
