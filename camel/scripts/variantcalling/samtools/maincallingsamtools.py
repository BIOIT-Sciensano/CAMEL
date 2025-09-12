#!/usr/bin/env python
import argparse
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional

from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger, initialize_logging
from camel.app.snakemake import snakemakeutils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import variant_calling


class MainCalling:
    """
    Class to run samtools variant calling using CAMEL.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        """
        self._args = MainCalling._parse_arguments(args)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        group = argument_parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--bam', type=Path, help="Input BAM file")
        argument_parser.add_argument('--reference', type=Path, required=True)
        argument_parser.add_argument('--reference-name')
        argument_parser.add_argument('--output', required=True)
        argument_parser.add_argument(
            '--output-consensus', help="If specified, the consensus sequence is saved in this file.")
        argument_parser.add_argument('--working-dir', type=Path, default=Path.cwd())
        argument_parser.add_argument('--ploidy', choices=['GRCh37', 'GRCh38', 'X', 'Y', '1'], default='1')
        argument_parser.add_argument('--calling-method', choices=('consensus', 'multiallelic'))
        argument_parser.add_argument('--skip-variants', choices=['snps', 'indels'])
        argument_parser.add_argument('--mutation-rate', type=float)
        argument_parser.add_argument('--minimal-mq', type=int)
        argument_parser.add_argument('--minimal-bq', type=int)
        argument_parser.add_argument('--output-all-sites', action='store_true')
        argument_parser.add_argument('--count-orphans', action='store_true')
        argument_parser.add_argument('--disable-baq', action='store_true')
        argument_parser.add_argument('--threads', type=int, default=8)
        argument_parser.add_argument('--input-type', choices=['illumina', 'ont'], default='illumina')
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the variant calling Snakefile to call the variants.
        :return: None
        """
        # Create configuration file
        config_data = self.__create_snakemake_config_data()
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

        # Copy input BAM file to the right location
        target_dir = self._args.working_dir / 'variant_calling' / 'read_mapping' / self._args.input_type
        if not target_dir.exists():
            target_dir.mkdir(parents=True)
        snakemakeutils.dump_object([ToolIOFile(self._args.bam)], target_dir / 'bam.io')

        # Run Snakemake to generate output file
        SnakePipelineUtils.run_snakemake(
            snakefile=variant_calling.SNAKEFILE,
            config_path=config_file,
            targets=[Path(variant_calling.OUTPUT_UNFILTERED_VCF)],
            working_dir=self._args.working_dir,
            threads=self._args.threads)

        # Generate consensus sequence
        if self._args.output_consensus:
            self.__generate_consensus_sequence(self._args.output_consensus, config_data)

        # Copy output
        logger.info("Collecting Snakemake output file")
        output_vcf_path = snakemakeutils.load_object(self._args.working_dir / variant_calling.OUTPUT_UNFILTERED_VCF)[0].path
        shutil.copyfile(output_vcf_path, self._args.output)

    def __create_snakemake_config_data(self) -> dict:
        """
        Creates a Snakemake configuration file.
        :return: Config file data
        """
        config_data = {
            'input_type': self._args.input_type,
            'reference': {
                'name': self._args.reference_name
                if self._args.reference_name
                else self._args.reference.name,
                'fasta': str(self._args.reference),
            },
            'sample_name': 'Sample',
            'variant_calling': {'ploidy': self._args.ploidy},
            'variant_filtering': {},
            'working_dir': str(self._args.working_dir),
        }
        for k in ['calling_method', 'skip_variants', 'mutation_rate', 'minimal_bq', 'minimal_mq', 'count_orphans',
                  'disable_baq']:
            if (k in self._args) and (vars(self._args)[k] is not None):
                config_data['variant_calling'][k] = vars(self._args)[k]
        if ('output_all_sites' in self._args) and self._args.output_all_sites is True:
            config_data['variant_calling']['variants_only'] = False
        return config_data

    def __generate_consensus_sequence(self, output_path: Path, config_data: dict[str, Any]) -> None:
        """
        Generates the consensus sequence by applying the detected variants to the reference sequence.
        :param output_path: Output path to save the consensus sequence
        :param config_data: Snakemake config data
        :return: None
        """
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir, 'consensus.yml')
        SnakePipelineUtils.run_snakemake(
            variant_calling.SNAKEFILE, config_file, [Path(variant_calling.OUTPUT_CONSENSUS)],  self._args.working_dir,
            self._args.threads)
        fasta_consensus = snakemakeutils.load_object(self._args.working_dir / variant_calling.OUTPUT_CONSENSUS)[0].path
        shutil.copyfile(fasta_consensus, output_path)


def main(args: Sequence[str] | None = None) -> None:
    """
    Entry point for the common interface.
    :param args: Command line arguments
    :return: None
    """
    script = MainCalling(args)
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
