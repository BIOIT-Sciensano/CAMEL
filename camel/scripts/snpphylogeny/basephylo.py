import abc
import dataclasses
import sys
from pathlib import Path

import click

from camel.app.core.errors import ToolExecutionError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.model import BaseOptions
from camel.app.toolkits.phylogeny import snpphylogenyutils
from camel.app.toolkits.phylogeny.snpphylogenyutils import (
    MappingInput,
    Sample, PhyloInput, PhyloOutput,
)
from camel.app.tools.mega.modelselection import ModelSelection


@dataclasses.dataclass(frozen=True)
class CommonOptions(BaseOptions):
    """
    Common options for the phylogeny pipelines.
    """
    threads: int = dataclasses.field(default=1, metadata={'help': 'Number of threads', 'show_default': True})
    working_dir: Path = dataclasses.field(default=Path.cwd(), metadata={'help': 'Working directory'})
    trim_reads: bool = dataclasses.field(default=False)
    missing_data: str = dataclasses.field(default='use_all_sites', metadata={
        'choices': ['complete_deletion', 'use_all_sites', 'partial_deletion'],
        'show_default': True,
    })
    site_cov_cutoff: int = dataclasses.field(default=10, metadata={
        'help': 'Minimum site coverage',
        'type': click.IntRange(0, 100)
    })
    bootstraps: int = dataclasses.field(default=100, metadata={
        'help': 'Number of bootstraps',
        'show_default': True,
    })
    ml_method: str = dataclasses.field(default='spr3', metadata={
        'choices': ['nni', 'spr3', 'spr5'],
        'show_default': True,
    })
    branch_swap: str = dataclasses.field(default='moderate', metadata={
        'choices': ['none', 'weak', 'very_weak', 'moderate', 'strong', 'very_strong'],
        'show_default': True,
    })
    include_ref: bool = dataclasses.field(default=False)
    report_include_bam: bool = dataclasses.field(default=False)


class BasePhylo(BaseScript[PhyloInput, PhyloOutput, CommonOptions], metaclass=abc.ABCMeta):
    """
    Base class for the SNP phylogeny pipelines.
    """

    def __init__(
            self, pipeline_name: str, version: str, script_in: PhyloInput, script_out: PhyloOutput, opts: CommonOptions
        ) -> None:
        """
        Initializes the base class.
        """
        self._pipeline_name = pipeline_name
        super().__init__(
            name=pipeline_name,
            version=version,
            script_in=script_in,
            script_out=script_out,
            script_opts=opts
        )
        self._report = snpphylogenyutils.initialize_report(self._pipeline_name, script_in, script_out)
        self._informs: list[dict] = []

    @property
    def samples_by_name(self) -> dict[str, Sample]:
        """
        Returns the samples as a dictionary with the sample name as key.
        :return: Samples by name
        """
        return {s.name_valid: s for s in self._script_in.samples}

    def _get_mapping_input(self) -> dict[Sample, MappingInput]:
        """
        Returns the input for the read mapping.
        :return: Mapping input per sample
        """
        logger.info("Preparing input for read mapping")
        fq_by_sample = snpphylogenyutils.symlink_input_files(self._script_in.samples, self._script_opts.working_dir)
        if self._script_opts.trim_reads:
            trimming_output_by_sample = snpphylogenyutils.trim_all_reads(
                fq_by_sample, self._script_opts.working_dir / 'trimming', self._script_opts.threads)
            snpphylogenyutils.add_trimming_section(self._report, trimming_output_by_sample)
            mapping_input_by_sample = {}
            for sample, output in trimming_output_by_sample.items():
                mapping_input_by_sample[sample] = MappingInput(
                    pe=output.trimmed_reads_pe,
                    se_fwd=output.trimmed_reads_se_fwd[0] if len(output.trimmed_reads_se_fwd) > 0 else None,
                    se_rev=output.trimmed_reads_se_rev[0] if len(output.trimmed_reads_se_rev) > 0 else None
                )
            self._informs.append(trimming_output_by_sample[self._script_in.samples[0]].informs_trimming)
            return mapping_input_by_sample
        else:
            snpphylogenyutils.add_trimming_section_empty(self._report)
            return {s: MappingInput(pe=fq) for s, fq in fq_by_sample.items()}

    def _run_model_selection(self, snp_matrix: ToolIOFile) -> ModelSelection:
        """
        Runs the model selection.
        Quits the programs if the SNP matrix is too small.
        :param snp_matrix: SNP matrix
        :return: Model selection tool instance
        """
        try:
            snpphylogenyutils.check_snp_matrix_size(snp_matrix.path)
        except ValueError as err:
            snpphylogenyutils.add_model_selection_section(self._report, error_message=str(err))
            sys.exit(0)
        else:
            model_selection = snpphylogenyutils.run_model_selection(
                snp_matrix=snp_matrix,
                dir_=self._script_opts.working_dir,
                missing_data=self._script_opts.missing_data,
                branch_swap=self._script_opts.branch_swap,
                site_cov_cutoff=self._script_opts.site_cov_cutoff,
                threads=self._script_opts.threads,
            )
            snpphylogenyutils.add_model_selection_section(self._report, model_selection=model_selection)
            self._informs.append(model_selection.informs)
            return model_selection

    def _run_tree_building(self, snp_matrix: ToolIOFile, model_selection: ModelSelection) -> None:
        """
        Runs the tree building.
        :param snp_matrix: SNP matrix
        :param model_selection: Model selection instance.
        :return: None
        """
        try:
            tree_building = snpphylogenyutils.run_tree_building(
                snp_matrix=snp_matrix,
                dir_=self._script_opts.working_dir,
                model=model_selection.informs['model'],
                rates=model_selection.informs['rates_among_sites'],
                bootstraps=self._script_opts.bootstraps,
                ml_method=self._script_opts.ml_method,
                branch_swap=self._script_opts.branch_swap,
                missing_data=self._script_opts.missing_data,
                site_cov_cutoff=self._script_opts.site_cov_cutoff,
                threads=self._script_opts.threads,
            )
            snpphylogenyutils.add_tree_building_section(self._report, tree_building.tool_outputs['NWK'][0].path)
            self._informs.append(tree_building.informs)
        except ToolExecutionError:
            snpphylogenyutils.add_tree_building_section(
                self._report, error_message='Error constructing tree, SNP matrix might be too small')
