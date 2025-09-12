#!/usr/bin/env python
import argparse
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional

from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.phylogeny.snpphylogenyutils import (
    MappingInput,
    Sample,
    SnpPhylogenyUtils,
)
from camel.app.components.workflows.variantcallingwrapper import VariantCallingOutput
from camel.app.components.workflows.variantfilteringwrapper import (
    VariantFilteringOutput,
)
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake import snakemakeutils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.bowtie2.bowtie2index import Bowtie2Index
from camel.scripts.snpphylogeny import (
    SNAKEFILE_SAMTOOLS_CALLING_ALL,
    SNAKEFILE_SAMTOOLS_FILTERING_ALL,
)
from camel.scripts.snpphylogeny.basephylo import BasePhylo
from camel.scripts.snpphylogeny.snakefile.samtools_calling_all import OUTPUT_CALLING_ALL
from camel.scripts.snpphylogeny.snakefile.samtools_filtering_all import (
    OUTPUT_FILTERING_ALL,
)


class MainSamtoolsPhylo(BasePhylo):
    """
    This class is used as the main script for the samtools SNP phylogeny.
    """

    # Type aliases
    FilteringOutBySample = dict[Sample, VariantFilteringOutput]
    CallingOutBySample = dict[Sample, VariantCallingOutput]

    PARAMETER_MAPPING = {
        'depth': {
            'min_total_depth': {'title': 'Min total depth', 'arg': lambda args: args.min_total_depth},
            'min_fwd_depth': {'title': 'Min forward depth', 'arg': lambda args: args.min_forward_depth},
            'min_rev_depth': {'title': 'Min reverse depth', 'arg': lambda args: args.min_reverse_depth}},
        'snp_quality': {
            'min_snp_quality': {'title': 'Min SNP quality', 'arg': lambda args: args.min_snp_quality}},
        'mapping_quality': {
            'min_mapping_quality': {'title': 'Min mapping quality', 'arg': lambda args: args.min_mapping_quality}},
        'distance': {
            'min_distance': {'title': 'Min variant distance', 'arg': lambda args: args.min_distance},
            'keep_best': {'title': 'Keep best', 'arg': lambda args: args.keep_best}},
        'zscore': {
            'min_zscore': {'title': 'Min Z-score', 'arg': lambda args: args.min_zscore},
            'y_multiplier': {'title': 'Y-multiplier', 'arg': lambda args: args.y_mult}}
    }

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        :param args: Main script arguments (optional))
        """
        super().__init__('Samtools', args)

    def run(self) -> None:
        """
        Runs the SNP phylogeny workflow.
        :return: None
        """
        # Prepare input
        mapping_input = self._get_mapping_input()
        reference = self.__prepare_reference()

        # Do variant calling and filtering
        calling_out_by_sample = self.__run_variant_calling_workflow(reference, mapping_input)
        filtering_out_by_sample = self.__run_variant_filtering_workflow(calling_out_by_sample)

        # Create SNP matrix
        snp_matrix = SnpPhylogenyUtils.construct_snp_matrix(
            [s.name_valid for s in self._samples],
            [filtering_out_by_sample[s].vcf_filtered for s in self._samples],
            self._args.working_dir / 'snp_matrix',
            self._args.include_ref,
            self._args.soft_filter
        )
        if self._args.output_fasta is not None:
            shutil.copyfile(snp_matrix.path, self._args.output_fasta)

        # Add sections to the report
        self.__add_snp_filtering_section(filtering_out_by_sample)
        self.__add_metrics_section(calling_out_by_sample, filtering_out_by_sample)
        self.__add_output_files_section(snp_matrix.path, calling_out_by_sample, filtering_out_by_sample)

        # Perform model selection & tree building
        model_selection = self._run_model_selection(snp_matrix)
        self._run_tree_building(snp_matrix, model_selection)

        # Add commands section
        all_informs = self._informs + calling_out_by_sample[self._samples[0]].informs_all + filtering_out_by_sample[
            self._samples[0]].informs
        self._report.add_html_object(SnakePipelineUtils.create_commands_section(all_informs, self._args.working_dir))
        self._report.save()

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        SnpPhylogenyUtils.add_common_arguments(argument_parser)
        argument_parser.add_argument('--ploidy', default='haploid', choices=['haploid', 'diploid'])
        argument_parser.add_argument('--calling-method', default='consensus', choices=['consensus', 'multiallelic'])
        argument_parser.add_argument('--soft-filter', action='store_true')
        argument_parser.add_argument('--min-total-depth', default=10, type=int)
        argument_parser.add_argument('--min-forward-depth', default=1, type=int)
        argument_parser.add_argument('--min-reverse-depth', default=1, type=int)
        argument_parser.add_argument('--min-snp-quality', default=25, type=int)
        argument_parser.add_argument('--min-mapping-quality', default=30, type=int)
        argument_parser.add_argument('--min-distance', default=10, type=int)
        argument_parser.add_argument('--keep-best', action='store_true')
        argument_parser.add_argument('--min-zscore', default=1.96, type=float)
        argument_parser.add_argument('--y-mult', default=10, type=float)
        argument_parser.add_argument('--include-ref', action='store_true')
        argument_parser.add_argument('--report-include-bam', action='store_true')
        return argument_parser.parse_args(args)

    def __prepare_reference(self) -> Path:
        """
        Prepares the reference genome.
        :return: Indexed reference genome prefix
        """
        dir_ref = self._args.working_dir / 'ref'
        if not dir_ref.is_dir():
            dir_ref.mkdir(parents=True)
        link_path = dir_ref / FileSystemHelper.make_valid(
            self._args.reference_name if self._args.reference_name else self._args.reference.name)
        if link_path.is_symlink():
            link_path.unlink()
        link_path.symlink_to(self._args.reference)
        bt2_index = Bowtie2Index()
        bt2_index.add_input_files({'FASTA_REF': [ToolIOFile(link_path)]})
        bt2_index.run(dir_ref)
        return bt2_index.tool_outputs['INDEX_GENOME_PREFIX'][0].value

    def __run_variant_calling_workflow(self, reference: Path, mapping_input: dict[Sample, MappingInput]) -> \
            CallingOutBySample:
        """
        Runs the variant filtering workflow in parallel on all samples.
        :param reference: Reference genome path
        :param mapping_input: Mapping input by sample
        :return: Dictionary with the output for each sample
        """
        working_dir = self._args.working_dir / 'calling'
        config_data = {
            'working_dir': str(working_dir),
            'samples': {s.name_valid: v.as_dict() for s, v in mapping_input.items()},
            'reference_info': {'name': self._args.reference_name, 'path': str(reference)},
            'options': {
                'ploidy': 1 if self._args.ploidy == 'haploid' else False,
                'calling_method': self._args.calling_method,
                'skip_variants': 'indels'
            }
        }
        config_file = SnakePipelineUtils.generate_config_file(config_data, working_dir)
        SnakePipelineUtils.run_snakemake(
            SNAKEFILE_SAMTOOLS_CALLING_ALL, config_file, [Path(OUTPUT_CALLING_ALL)], working_dir, self._args.threads)
        return {self.samples_by_name[name]: output for name, output in snakemakeutils.load_object(working_dir / OUTPUT_CALLING_ALL).items()}

    def __run_variant_filtering_workflow(self, calling_output_by_sample: CallingOutBySample) -> FilteringOutBySample:
        """
        Runs the variant filtering workflow.
        :return: Dictionary with the output for each sample
        """
        working_dir = self._args.working_dir / 'filtering'
        samples = {
            s.name_valid: {'VCF': str(o.vcf_unfiltered.path), 'BAM': str(o.bam_file.path)} for s, o in
            calling_output_by_sample.items()
        }
        config_data = {'working_dir': str(working_dir), 'samples': samples, 'options': self.__get_filtering_options()}
        config_file = SnakePipelineUtils.generate_config_file(config_data, working_dir)
        output_path = OUTPUT_FILTERING_ALL
        SnakePipelineUtils.run_snakemake(
            SNAKEFILE_SAMTOOLS_FILTERING_ALL, config_file, [output_path], working_dir, self._args.threads)
        return {self.samples_by_name[name]: output for name, output in snakemakeutils.load_object(working_dir / output_path).items()}

    def __get_filtering_options(self) -> dict[str, Any]:
        """
        Returns the dictionary with filtering options.
        :return: Filtering options as a dictionary
        """
        options = {}
        if self._args.soft_filter:
            options['soft_filter'] = True
        for group, params in MainSamtoolsPhylo.PARAMETER_MAPPING.items():
            options[group] = {name: value['arg'](self._args) for name, value in params.items()}
        return options

    def __add_snp_filtering_section(self, filtering_out_by_sample: FilteringOutBySample) -> None:
        """
        Adds the section with the SNP filtering results to the report.
        :return: None
        """
        section = HtmlReportSection('SNP filtering')
        section.add_paragraph('This table show the number of SNPs that passed each of the filtering steps.')
        table_data = []
        for sample, output in filtering_out_by_sample.items():
            stats = output.stats
            row = [sample.name_full]
            for f in ('depth', 'snp_qual', 'mapping_qual', 'distance', 'zscore'):
                row.append('{}/{}'.format(stats[f]['variants_out'], stats[f]['variants_in']))
            table_data.append(row)
        header = ['Sample', 'Depth', 'SNP quality', 'Mapping quality', 'Distance', 'Z-score']
        section.add_table(table_data, header, [('class', 'data')])
        self.__add_filtering_options_table(section)
        if self._args.soft_filter:
            section.add_alert(
                "Soft filtering enabled: filtered variants are kept in the VCF file and annotated in the FILTER column",
                'info')
        self._report.add_html_object(section)
        self._report.save()

    def __add_filtering_options_table(self, section: HtmlReportSection) -> None:
        """
        Adds a table with the filtering options to the report.
        :return: None
        """
        section.add_paragraph('Filtering settings:')
        table_data = []
        for group, params in MainSamtoolsPhylo.PARAMETER_MAPPING.items():
            for param_name, param_value in params.items():
                table_data.append([param_value['title'], param_value['arg'](self._args)])
        header = ['Filter', 'Value']
        section.add_table(table_data, header, [('class', 'data')])

    def __add_metrics_section(self, calling_out_by_sample: CallingOutBySample,
                              filtering_out_by_sample: FilteringOutBySample) -> None:
        """
        Adds the section with the analysis metrics.
        :return: None
        """
        stats = {sample: [
            2 * int(calling_out_by_sample[sample].informs_mapping['stats_paired_reads_in'].split(' ')[0]) +
            int(calling_out_by_sample[sample].informs_mapping.get('stats_singe_reads_in', "0").split(' ')[0]),
            calling_out_by_sample[sample].informs_mapping['stats_map_rate'],
            calling_out_by_sample[sample].nb_of_variants,
            filtering_out_by_sample[sample].nb_of_variants] for sample in self._samples}
        header = ['Sample', 'Total reads', 'Mapping rate (%)', 'Nb. of SNPs (unfiltered)', 'Nb. of SNPs (filtered)']
        SnpPhylogenyUtils.add_metrics_section(self._report, stats, header)

    def __add_output_files_section(self, snp_matrix: Path, calling_out_by_sample: CallingOutBySample,
                                   filtering_out_by_sample: FilteringOutBySample) -> None:
        """
        Adds the section with the output files.
        :param snp_matrix: SNP matrix path
        :return: None
        """
        output_files = {sample: [
            calling_out_by_sample[sample].bam_file.path if self._args.report_include_bam else None,
            calling_out_by_sample[sample].vcf_unfiltered.path,
            filtering_out_by_sample[sample].vcf_filtered.path
        ] for sample in self._samples}
        column_names = ['Alignment (BAM)', 'SNPs unfiltered (VCF)', 'SNPs filtered (VCF)']
        SnpPhylogenyUtils.add_output_files_section(self._report, column_names, output_files, snp_matrix)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainSamtoolsPhylo()
    main.run()
