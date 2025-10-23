#!/usr/bin/env python
import argparse
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional

from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.utils import fileutils
from camel.app.loggers import logger, initialize_logging
from camel.app.tools.medaka.medakainference import MedakaInference
from camel.app.tools.medaka.medakavcf import MedakaVcf
from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex


class MainCallingMedaka:
    """
    This class contains the main script for the Medaka variant calling tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        """
        self._args = MainCallingMedaka._parse_arguments(args)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--bam', type=Path, help="Input BAM file", required=True)
        argument_parser.add_argument('--reference', help="The reference fasta file to use", type=Path, required=True)
        argument_parser.add_argument('--reference-name')
        argument_parser.add_argument('--output', help="Define custom output file in stead of default generated one")
        argument_parser.add_argument('--working-dir', help='Working directory', type=Path, default=Path.cwd())
        argument_parser.add_argument('--threads', type=int, default=4)
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the tool.
        :return: None
        """
        input_dict = self.__prepare_input()

        # Index BAM file
        logger.info('Using samtools to index input BAM file')
        samtools_index = SamtoolsIndex()
        samtools_index.add_input_files(input_dict)
        samtools_index.run(input_dict['BAM'][0].path.parent)

        # Run Medaka Inference
        medaka_inference = MedakaInference()
        logger.info(f'Running {medaka_inference.name} to get hdf file')
        medaka_inference.add_input_files(input_dict)
        medaka_inference.update_parameters(threads=self._args.threads)
        hdf_base = Path(self._args.bam).stem
        hdf_file = self._args.working_dir / Path(hdf_base + '.hdf')
        medaka_inference.update_parameters(output=hdf_file)
        medaka_inference.run(self._args.working_dir)

        # Run Medaka VCF
        medaka_vcf = MedakaVcf()
        logger.info(f'Running {medaka_vcf.name} to call the variants based on reference fasta and hdf file')
        medaka_vcf.add_input_files({'HDF': [ToolIOFile(Path(hdf_file))], 'FASTA': input_dict['FASTA']})
        if self._args.output is not None:
            medaka_vcf.update_parameters(output=self._args.output)
        else:
            vcf_base = Path(hdf_file).stem
            vcf_file = self._args.working_dir / Path(vcf_base + '.vcf')
            medaka_vcf.update_parameters(output=vcf_file)
        medaka_vcf.run(self._args.working_dir)

    def __prepare_input(self) -> dict[str, Any]:
        """
        Prepares the input for the Medaka variant tool.
        :return: Input dictionary
        """
        self._args.working_dir.mkdir(parents=True, exist_ok=True)
        input_dict: dict[str, None | list[ToolIOFile]] = {'BAM': None, 'FASTA': None}

        # input BAM file
        bam_file = Path(self._args.bam)
        path_bam = self._args.working_dir / fileutils.make_valid(bam_file.name)
        path_bam.symlink_to(bam_file)
        input_dict['BAM'] = [ToolIOFile(path_bam)]

        # Reference genome
        reference_file = Path(self._args.reference)
        path_ref = self._args.working_dir / fileutils.make_valid(reference_file.name)
        path_ref.symlink_to(reference_file)
        input_dict['FASTA'] = [ToolIOFile(path_ref)]

        return input_dict


if __name__ == '__main__':
    initialize_logging()
    main = MainCallingMedaka()
    main.run()
