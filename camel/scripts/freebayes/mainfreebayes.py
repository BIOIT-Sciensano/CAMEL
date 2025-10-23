#!/usr/bin/env python
import argparse
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import initialize_logging
from camel.app.tools.freebayes.freebayes import Freebayes


class MainFreebayesCalling:
    """
    Class to run freebayes variant calling using CAMEL.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        """
        self._args = MainFreebayesCalling._parse_arguments(args)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--bam', type=Path, help="Input BAM file", required=True)
        argument_parser.add_argument('--reference', type=Path, required=True)
        argument_parser.add_argument('--reference-name')
        argument_parser.add_argument('--output', required=True)
        argument_parser.add_argument('--working-dir', type=Path, default=Path.cwd())
        argument_parser.add_argument('--min-base-quality', type=int, default=0)
        argument_parser.add_argument('--min-coverage', type=int, default=0)
        argument_parser.add_argument('--min-mapping-quality', type=int, default=1)
        argument_parser.add_argument('--min-supporting-allele-qsum', type=int, default=0)
        argument_parser.add_argument('--ploidy', required=True, type=int, default=1)
        argument_parser.add_argument('--report-monomorphic', action="store_true", default=False)
        argument_parser.add_argument('--standard-filters', action="store_true", default=False)
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        self.__run_freebayes()

    def __run_freebayes(self) -> Freebayes:
        """
        Runs Freebayes.
        :return: Freebayes tool instance.
        """
        freebayes = Freebayes()
        if self._args.reference is not None:
            freebayes.add_input_files({'FASTA': [ToolIOFile(self._args.reference)]})
        if self._args.bam is not None:
            freebayes.add_input_files({'BAM': [ToolIOFile(self._args.bam)]})

        if self._args.standard_filters:
            freebayes.update_parameters(standard_filters=True)
        else:
            if self._args.min_base_quality:
                freebayes.update_parameters(min_base_quality=self._args.min_base_quality)
            if self._args.min_mapping_quality != 1:
                freebayes.update_parameters(min_mapping_quality=self._args.min_mapping_quality)
            if self._args.min_supporting_allele_qsum:
                freebayes.update_parameters(min_supporting_allele_qsum=self._args.min_supporting_allele_qsum)

        if self._args.report_monomorphic:
            freebayes.update_parameters(report_monomorphic=True)
        if self._args.min_coverage:
            freebayes.update_parameters(min_coverage=self._args.min_coverage)

        freebayes.update_parameters(ploidy=self._args.ploidy, vcf=self._args.output)

        freebayes.run(self._args.working_dir)
        return freebayes


if __name__ == '__main__':
    initialize_logging()
    freebayes_main = MainFreebayesCalling()
    freebayes_main.run()
