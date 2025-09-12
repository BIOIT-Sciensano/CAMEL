#!/usr/bin/env python
import argparse
import itertools
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

from camel.app.camel import Camel
from camel.app.components.phylogeny.snpphylogenyutils import SnpPhylogenyUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.cfsan.cfsansnppipeline import CfsanSnpPipeline
from camel.scripts.snpphylogeny.basephylo import BasePhylo


class MainCfsanPhylo(BasePhylo):
    """
    This class is used as the main script for the CFSAN SNP phylogeny.
    """

    def __init__(self, args: Optional[Sequence[str]] = None):
        """
        Initializes the main script.
        :param args: Main script arguments (optional))
        """
        super().__init__('CFSAN', args)

    def run(self) -> None:
        """
        Runs the SNP phylogeny workflow.
        :return: None
        """
        # Prepare input
        mapping_input = self._get_mapping_input()
        reference = self.__prepare_reference()

        # Run SNP calling
        fq_files = list(itertools.chain.from_iterable([mi.reads_raw for mi in mapping_input]))
        cfsan = self.__run_cfsan(reference, fq_files)

        # Get SNP matrix
        snp_matrix = cfsan.tool_outputs['FASTA_Preserved' if self._args.selected_matrix == 'preserved' else 'FASTA'][0]
        if self._args.output_fasta is not None:
            shutil.copyfile(snp_matrix.path, self._args.output_fasta)

        # Add sections to the report
        self.__add_metrics_section(cfsan)
        self.__add_output_files_section(cfsan, snp_matrix)

        # Perform model selection & tree building
        model_selection = self._run_model_selection(snp_matrix)
        self._run_tree_building(snp_matrix, model_selection)

        # Add commands section
        all_informs = self._informs
        self._report.add_html_object(SnakePipelineUtils.create_commands_section(all_informs, self._args.working_dir))
        self._report.save()

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        SnpPhylogenyUtils.add_common_arguments(argument_parser)
        argument_parser.add_argument('--selected-matrix', choices=['preserved', 'regular'])
        argument_parser.add_argument('--report-include-bam', action='store_true')
        return argument_parser.parse_args(args)

    def __prepare_reference(self) -> Path:
        """
        Prepares the reference genome.
        :return: Path to reference genome.
        """
        reference_name = self._args.reference_name if self._args.reference_name is not None else \
            self._args.reference.name
        reference_path = self._args.working_dir / reference_name
        if not reference_path.exists():
            reference_path.symlink_to(self._args.reference)
        return reference_path

    def __run_cfsan(self, reference: Path, mapping_input: list[ToolIOFile]) -> CfsanSnpPipeline:
        """
        Runs the CFSAN SNP pipeline.
        :param reference: Reference genome
        :param mapping_input: Mapping input dictionary
        :return: CFSAN tool instance
        """
        cfsan_input = {
            'VAL_name': [ToolIOValue(s.name_valid) for s in self._samples],
            'FASTQ': mapping_input,
            'FASTA': [ToolIOFile(reference)]
        }
        cfsan = CfsanSnpPipeline()
        cfsan.add_input_files(cfsan_input)
        cfsan.run(self._args.working_dir)
        self._informs.append(cfsan.informs)
        return cfsan

    def __add_metrics_section(self, cfsan: CfsanSnpPipeline) -> None:
        """
        Adds the metrics to the output report.
        :param cfsan: CFSAN instance
        :return: None
        """
        stats = {s: [
            cfsan.informs[s.name_valid]['Number_of_Reads'],
            cfsan.informs[s.name_valid]['Percent_of_Reads_Mapped'],
            cfsan.informs[s.name_valid]['Average_Pileup_Depth'],
            cfsan.informs[s.name_valid]['Phase1_SNPs'],
            cfsan.informs[s.name_valid]['Phase1_Preserved_SNPs']
        ] for s in self._samples}
        header = ['Sample', 'Total reads', 'Mapping rate (%)', 'Avg. pileup depth', 'Nb. of SNPs (unfiltered)',
                  'Nb. of SNPs (filtered)']
        SnpPhylogenyUtils.add_metrics_section(self._report, stats, header)

    def __add_output_files_section(self, cfsan: CfsanSnpPipeline, snp_matrix: ToolIOFile) -> None:
        """
        Adds the section with the output files to the report.
        :param cfsan: CFSAN instance
        :param snp_matrix: SNP matrix
        :return: None
        """
        output_files = {self._samples[i]: [
            cfsan.tool_outputs['BAM'][i].path if self._args.report_include_bam else None,
            cfsan.tool_outputs['VCF'][i].path,
            cfsan.tool_outputs['VCF_preserved'][i].path
        ] for i in range(0, len(cfsan.tool_outputs['Sample_names']))}
        header = ['Alignment (BAM)', 'SNPs filtered (VCF)', 'SNPs filtered preserved (VCF)']
        SnpPhylogenyUtils.add_output_files_section(self._report, header, output_files, snp_matrix.path)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainCfsanPhylo()
    main.run()
