#!/usr/bin/env python
import dataclasses
import itertools
import shutil
from pathlib import Path

import click

from camel.app.cli import cliutils
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.reports import reportutils
from camel.app.loggers import initialize_logging
from camel.app.scriptutils.model import BaseOptions
from camel.app.toolkits.phylogeny import snpphylogenyutils
from camel.app.toolkits.phylogeny.snpphylogenyutils import PhyloInput, PhyloOutput
from camel.app.tools.cfsan.cfsansnppipeline import CfsanSnpPipeline
from camel.scripts.snpphylogeny.basephylo import BasePhylo, CommonOptions


@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Specific options for the CFSAN SNP phylogeny pipeline.
    """
    selected_matrix: str = dataclasses.field(default='preserved', metadata={'choices': ['preserved', 'regular']})

class MainCfsanPhylo(BasePhylo):
    """
    This class is used as the main script for the CFSAN SNP phylogeny.
    """

    def __init__(self, in_: PhyloInput, out_: PhyloOutput, opts_common: CommonOptions, opts_custom: Options) -> None:
        """
        Initializes the main script.
        """
        super().__init__(
            pipeline_name='SNP phylogeny (CFSAN)',
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

        # Run SNP calling
        fq_files = list(itertools.chain.from_iterable([mi.reads_raw for mi in mapping_input]))
        cfsan = self.__run_cfsan(reference, fq_files)

        # Get SNP matrix
        snp_matrix = cfsan.tool_outputs['FASTA_Preserved' if self._opts_custom.selected_matrix == 'preserved' else 'FASTA'][0]
        if self._script_out.output_fasta is not None:
            shutil.copyfile(snp_matrix.path, self._script_out.output_fasta)

        # Add sections to the report
        self.__add_metrics_section(cfsan)
        self.__add_output_files_section(cfsan, snp_matrix)

        # Perform model selection and tree building
        model_selection = self._run_model_selection(snp_matrix)
        self._run_tree_building(snp_matrix, model_selection)

        # Add the commands section
        all_informs = self._informs
        self._report.add_html_object(reportutils.create_commands_section(all_informs, self._script_opts.working_dir))
        self._report.save()

    def __prepare_reference(self) -> Path:
        """
        Prepares the reference genome.
        :return: Path to reference genome.
        """
        reference_name = self._script_in.reference_name if self._script_in.reference_name is not None else \
            self._script_in.reference.name
        reference_path = self._script_opts.working_dir / reference_name
        if not reference_path.exists():
            reference_path.symlink_to(self._script_in.reference)
        return reference_path

    def __run_cfsan(self, reference: Path, mapping_input: list[ToolIOFile]) -> CfsanSnpPipeline:
        """
        Runs the CFSAN SNP pipeline.
        :param reference: Reference genome
        :param mapping_input: Mapping input dictionary
        :return: CFSAN tool instance
        """
        cfsan_input = {
            'VAL_name': [ToolIOValue(s.name_valid) for s in self._script_in.samples],
            'FASTQ': mapping_input,
            'FASTA': [ToolIOFile(reference)]
        }
        cfsan = CfsanSnpPipeline()
        cfsan.add_input_files(cfsan_input)
        cfsan.run(self._script_opts.working_dir)
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
        ] for s in self._script_in.samples}
        header = ['Sample', 'Total reads', 'Mapping rate (%)', 'Avg. pileup depth', 'Nb. of SNPs (unfiltered)',
                  'Nb. of SNPs (filtered)']
        snpphylogenyutils.add_metrics_section(self._report, stats, header)

    def __add_output_files_section(self, cfsan: CfsanSnpPipeline, snp_matrix: ToolIOFile) -> None:
        """
        Adds the section with the output files to the report.
        :param cfsan: CFSAN instance
        :param snp_matrix: SNP matrix
        :return: None
        """
        output_files = {self._script_in.samples[i]: [
            cfsan.tool_outputs['BAM'][i].path if self._script_opts.report_include_bam else None,
            cfsan.tool_outputs['VCF'][i].path,
            cfsan.tool_outputs['VCF_preserved'][i].path
        ] for i in range(0, len(cfsan.tool_outputs['Sample_names']))}
        header = ['Alignment (BAM)', 'SNPs filtered (VCF)', 'SNPs filtered preserved (VCF)']
        snpphylogenyutils.add_output_files_section(self._report, header, output_files, snp_matrix.path)


@click.command(name='snp_phylogeny_cfsan', short_help='SNP phylogeny with CFSAN variant calling.')
@snpphylogenyutils.add_input_options
@cliutils.add_click_options_from_dataclass(PhyloOutput)
@cliutils.add_click_options_from_dataclass(CommonOptions)
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    SNP phylogeny pipeline with SAMtools variant calling.
    """
    script_in = snpphylogenyutils.parse_input_options(kwargs)
    script = MainCfsanPhylo(
        in_=script_in,
        out_=PhyloOutput(**cliutils.from_kwargs(PhyloOutput, kwargs)),
        opts_common=CommonOptions(**cliutils.from_kwargs(CommonOptions, kwargs)),
        opts_custom=Options(**cliutils.from_kwargs(Options, kwargs)),
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
