#!/usr/bin/env python
import dataclasses
import shutil
from pathlib import Path
from typing import Any

import click

from camel.app.cli import cliutils
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.reports import reportutils
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.snakemake import snakemakeutils
from camel.app.core.snakemake import snakepipelineutils
from camel.app.core.utils import fileutils
from camel.app.loggers import initialize_logging
from camel.app.scriptutils.model import BaseOptions
from camel.app.toolkits.phylogeny import snpphylogenyutils
from camel.app.toolkits.phylogeny.snpphylogenyutils import MappingInput, Sample, PhyloInput, PhyloOutput
from camel.app.tools.bowtie2.bowtie2index import Bowtie2Index
from camel.app.wrappers.variantcallingwrapper import VariantCallingOutput
from camel.app.wrappers.variantfilteringwrapper import VariantFilteringOutput
from camel.scripts.snpphylogeny import (
    SNAKEFILE_SAMTOOLS_CALLING_ALL,
    SNAKEFILE_SAMTOOLS_FILTERING_ALL,
)
from camel.scripts.snpphylogeny.basephylo import BasePhylo, CommonOptions
from camel.scripts.snpphylogeny.snakefile.samtools_calling_all import OUTPUT_CALLING_ALL
from camel.scripts.snpphylogeny.snakefile.samtools_filtering_all import (
    OUTPUT_FILTERING_ALL,
)


@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Specific options for the samtools SNP phylogeny pipeline.
    """
    ploidy: str = dataclasses.field(default='haploid', metadata={'choices': ['haploid', 'diploid']})
    calling_method: str = dataclasses.field(default='consensus', metadata={'choices': ['consensus', 'multiallelic']})
    soft_filter: bool = dataclasses.field(default=False)
    min_total_depth: int = dataclasses.field(default=10)
    min_fwd_depth: int = dataclasses.field(default=1)
    min_rev_depth: int = dataclasses.field(default=1)
    min_snp_quality: int = dataclasses.field(default=25)
    min_mapping_quality: int = dataclasses.field(default=30)
    min_distance: int = dataclasses.field(default=10)
    keep_best: bool = dataclasses.field(default=False)
    min_zscore: float = dataclasses.field(default=1.96)
    y_multiplier: float = dataclasses.field(default=10)


class MainSamtoolsPhylo(BasePhylo):
    """
    This class is used as the main script for the samtools SNP phylogeny.
    """

    # Type aliases
    FilteringOutBySample = dict[Sample, VariantFilteringOutput]
    CallingOutBySample = dict[Sample, VariantCallingOutput]

    PARAMETER_MAPPING = {
        "depth": {
            "min_total_depth": {"title": "Min total depth"},
            "min_fwd_depth": {"title": "Min forward depth"},
            "min_rev_depth": {"title": "Min reverse depth"},
        },
        "snp_quality": {
            "min_snp_quality": {"title": "Min SNP quality"}},
        "mapping_quality": {
            "min_mapping_quality": {"title": "Min mapping quality"}},
        "distance": {
            "min_distance": {"title": "Min variant distance"},
            "keep_best": {"title": "Keep best"},
        },
        "zscore": {
            "min_zscore": {"title": "Min Z-score"},
            "y_multiplier": {"title": "Y-multiplier"},
        },
    }

    def __init__(self, in_: PhyloInput, out_: PhyloOutput, opts_common: CommonOptions, opts_custom: Options) -> None:
        """
        Initializes the main script.
        """
        super().__init__(
            pipeline_name='SNP phylogeny (SAMtools)',
            version='1.0',
            script_in=in_,
            script_out=out_,
            opts=opts_common
        )
        self._opts_custom = opts_custom

    def _execute(self) -> None:
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
        snp_matrix = snpphylogenyutils.construct_snp_matrix(
            [s.name_valid for s in self._script_in.samples],
            [filtering_out_by_sample[s].vcf_filtered for s in self._script_in.samples],
            self._script_opts.working_dir / "snp_matrix",
            self._script_opts.include_ref,
            self._opts_custom.soft_filter,
        )
        if self._script_out.output_fasta is not None:
            shutil.copyfile(snp_matrix.path, self._script_out.output_fasta)

        # Add sections to the report
        self.__add_snp_filtering_section(filtering_out_by_sample)
        self.__add_metrics_section(calling_out_by_sample, filtering_out_by_sample)
        self.__add_output_files_section(
            snp_matrix.path, calling_out_by_sample, filtering_out_by_sample
        )

        # Perform model selection & tree building
        model_selection = self._run_model_selection(snp_matrix)
        self._run_tree_building(snp_matrix, model_selection)

        # Add commands section
        all_informs = (
            self._informs
            + calling_out_by_sample[self._script_in.samples[0]].informs_all
            + filtering_out_by_sample[self._script_in.samples[0]].informs
        )
        self._report.add_html_object(reportutils.create_commands_section(all_informs, self._script_opts.working_dir))
        self._report.save()

    def __prepare_reference(self) -> Path:
        """
        Prepares the reference genome.
        :return: Indexed reference genome prefix
        """
        dir_ref = self._script_opts.working_dir / "ref"
        if not dir_ref.is_dir():
            dir_ref.mkdir(parents=True)
        link_path = dir_ref / fileutils.make_valid(
            self._script_in.reference_name
            if self._script_in.reference_name
            else self._script_in.reference.name
        )
        if link_path.is_symlink():
            link_path.unlink()
        link_path.symlink_to(self._script_in.reference)
        bt2_index = Bowtie2Index()
        bt2_index.add_input_files({"FASTA_REF": [ToolIOFile(link_path)]})
        bt2_index.run(dir_ref)
        return bt2_index.tool_outputs["INDEX_GENOME_PREFIX"][0].value

    def __run_variant_calling_workflow(
        self, reference: Path, mapping_input: dict[Sample, MappingInput]
    ) -> CallingOutBySample:
        """
        Runs the variant filtering workflow in parallel on all samples.
        :param reference: Reference genome path
        :param mapping_input: Mapping input by sample
        :return: Dictionary with the output for each sample
        """
        working_dir = self._script_opts.working_dir / "calling"
        config_data = {
            "working_dir": str(working_dir),
            "samples": {s.name_valid: v.as_dict() for s, v in mapping_input.items()},
            "reference_info": {
                "name": self._script_in.reference_name,
                "path": str(reference),
            },
            "options": {
                "ploidy": 1 if self._opts_custom.ploidy == "haploid" else False,
                "calling_method": self._opts_custom.calling_method,
                "skip_variants": "indels",
            },
        }
        config_file = snakepipelineutils.generate_config_file(config_data, working_dir)
        snakepipelineutils.run_snakemake(
            snakefile=SNAKEFILE_SAMTOOLS_CALLING_ALL,
            config_path=config_file,
            targets=[Path(OUTPUT_CALLING_ALL)],
            working_dir=working_dir,
            threads=self._script_opts.threads,
        )
        return {
            self.samples_by_name[name]: output for name, output in snakemakeutils.load_object(
                working_dir / OUTPUT_CALLING_ALL).items()}

    def __run_variant_filtering_workflow(
        self, calling_output_by_sample: CallingOutBySample
    ) -> FilteringOutBySample:
        """
        Runs the variant filtering workflow.
        :return: Dictionary with the output for each sample
        """
        working_dir = self._script_opts.working_dir / "filtering"
        samples = {
            s.name_valid: {
                "VCF": str(o.vcf_unfiltered.path),
                "BAM": str(o.bam_file.path),
            }
            for s, o in calling_output_by_sample.items()
        }
        config_data = {
            "working_dir": str(working_dir),
            "samples": samples,
            "options": self.__get_filtering_options(),
        }
        config_file = snakepipelineutils.generate_config_file(config_data, working_dir)
        output_path = Path(OUTPUT_FILTERING_ALL)
        snakepipelineutils.run_snakemake(
            snakefile=SNAKEFILE_SAMTOOLS_FILTERING_ALL,
            config_path=config_file,
            targets=[output_path],
            working_dir=working_dir,
            threads=self._script_opts.threads,
        )
        return {
            self.samples_by_name[name]: output for name, output in snakemakeutils.load_object(
                working_dir / output_path).items()}

    def __get_filtering_options(self) -> dict[str, Any]:
        """
        Returns the dictionary with filtering options.
        :return: Filtering options as a dictionary
        """
        options = {}
        if self._opts_custom.soft_filter:
            options["soft_filter"] = True
        for group, params in MainSamtoolsPhylo.PARAMETER_MAPPING.items():
            options[group] = {name: getattr(self._opts_custom, name) for name, value in params.items()}
        return options

    def __add_snp_filtering_section(
        self, filtering_out_by_sample: FilteringOutBySample
    ) -> None:
        """
        Adds the section with the SNP filtering results to the report.
        :return: None
        """
        section = HtmlReportSection("SNP filtering")
        section.add_paragraph(
            "This table show the number of SNPs that passed each of the filtering steps."
        )
        table_data = []
        for sample, output in filtering_out_by_sample.items():
            stats = output.stats
            row = [sample.name_full]
            for f in ("depth", "snp_qual", "mapping_qual", "distance", "zscore"):
                row.append(f"{stats[f]['variants_out']}/{stats[f]['variants_in']}")
            table_data.append(row)
        header = [
            "Sample",
            "Depth",
            "SNP quality",
            "Mapping quality",
            "Distance",
            "Z-score",
        ]
        section.add_table(table_data, header, [("class", "data")])
        self.__add_filtering_options_table(section)
        if self._opts_custom.soft_filter:
            section.add_alert(
                "Soft filtering enabled: filtered variants are kept in the VCF file and annotated in the FILTER column",
                "info",
            )
        self._report.add_html_object(section)
        self._report.save()

    def __add_filtering_options_table(self, section: HtmlReportSection) -> None:
        """
        Adds a table with the filtering options to the report.
        :return: None
        """
        section.add_paragraph("Filtering settings:")
        table_data = []
        print(self._opts_custom)
        for group, params in MainSamtoolsPhylo.PARAMETER_MAPPING.items():
            for param_name, param_value in params.items():
                table_data.append([param_value["title"], getattr(self._opts_custom, param_name)])
        header = ["Filter", "Value"]
        section.add_table(table_data, header, [("class", "data")])

    def __add_metrics_section(
        self,
        calling_out_by_sample: CallingOutBySample,
        filtering_out_by_sample: FilteringOutBySample,
    ) -> None:
        """
        Adds the section with the analysis metrics.
        :return: None
        """
        stats = {
            sample: [
                2 * int(calling_out_by_sample[sample].informs_mapping["stats_paired_reads_in"].split(" ")[0])
                + int(
                    calling_out_by_sample[sample]
                    .informs_mapping.get("stats_singe_reads_in", "0")
                    .split(" ")[0]
                ),
                calling_out_by_sample[sample].informs_mapping["stats_map_rate"],
                calling_out_by_sample[sample].nb_of_variants,
                filtering_out_by_sample[sample].nb_of_variants,
            ]
            for sample in self._script_in.samples
        }
        header = [
            "Sample",
            "Total reads",
            "Mapping rate (%)",
            "Nb. of SNPs (unfiltered)",
            "Nb. of SNPs (filtered)",
        ]
        snpphylogenyutils.add_metrics_section(self._report, stats, header)

    def __add_output_files_section(
        self,
        snp_matrix: Path,
        calling_out_by_sample: CallingOutBySample,
        filtering_out_by_sample: FilteringOutBySample,
    ) -> None:
        """
        Adds the section with the output files.
        :param snp_matrix: SNP matrix path
        :return: None
        """
        output_files = {
            sample: [
                (
                    calling_out_by_sample[sample].bam_file.path
                    if self._script_opts.report_include_bam
                    else None
                ),
                calling_out_by_sample[sample].vcf_unfiltered.path,
                filtering_out_by_sample[sample].vcf_filtered.path,
            ]
            for sample in self._script_in.samples
        }
        column_names = [
            "Alignment (BAM)",
            "SNPs unfiltered (VCF)",
            "SNPs filtered (VCF)",
        ]
        snpphylogenyutils.add_output_files_section(
            self._report, column_names, output_files, snp_matrix
        )


@click.command(name='snp_phylogeny_samtools', short_help='SNP phylogeny with SAMtools variant calling.')
@snpphylogenyutils.add_input_options
@cliutils.add_click_options_from_dataclass(PhyloOutput)
@cliutils.add_click_options_from_dataclass(CommonOptions)
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    SNP phylogeny pipeline with SAMtools variant calling.
    """
    script_in = snpphylogenyutils.parse_input_options(kwargs)
    script = MainSamtoolsPhylo(
        in_=script_in,
        out_=PhyloOutput(**cliutils.from_kwargs(PhyloOutput, kwargs)),
        opts_common=CommonOptions(**cliutils.from_kwargs(CommonOptions, kwargs)),
        opts_custom=Options(**cliutils.from_kwargs(Options, kwargs)),
    )
    script.run()


if __name__ == "__main__":
    initialize_logging()
    main()
