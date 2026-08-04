#!/usr/bin/env python
import dataclasses
import shutil
from pathlib import Path
from typing import Any

import click
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.io.tooliovalue import ToolIOValue

from camel.app.cli import cliutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.loggers import initialize_logging
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.model import BaseInput, BaseOptions, BaseOutput
from camel.app.toolkits.phylogeny.megautils import MEGAUtils
from camel.app.tools.mega.mltreeconstruction import MLTreeConstruction
from camel.app.tools.mega.modelselection import ModelSelection
from camel.app.tools.snpmatrix.snpmatrixconstructor import SnpMatrixConstructor
from camel.version import __VERSION__


@dataclasses.dataclass(frozen=True)
class Input(BaseInput):
    """
    Input data for the MEGA script.
    """
    fasta: Path | None = None
    vcf: list[tuple[Path, str]] = dataclasses.field(default_factory=list)

    def validate(self) -> bool:
        """
        Checks if the provided input is valid.
        :return: True if valid, False otherwise
        """
        if self.fasta is None and len(self.vcf) == 0:
            raise InvalidToolInputError('Either FASTA or VCF must be provided')
        return True


@dataclasses.dataclass(frozen=True)
class Output(BaseOutput):
    """
    Output data for the MEGA script.
    """
    output_tree: Path | None = None
    output_model: Path | None = None
    output_snp_matrix: Path | None = None


@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Options for the MEGA script.
    """
    working_dir: Path = dataclasses.field(default=Path.cwd(), metadata={'help': 'Working directory'})
    action: str | None = dataclasses.field(default='both', metadata={
        'choices': ['both', 'tree', 'model'], 'show_default': True})
    model: str | None = dataclasses.field(default='JC', metadata={
        'choices': ['JC', 'K2', 'T92', 'HKY', 'TN93', 'GTR'], 'show_default': True})
    missing_data: str | None = dataclasses.field(default='use_all_sites', metadata={
            'choices': ['use_all_sites', 'complete_deletion', 'partial_deletion'],
            'show_default': True
        })
    branch_swap: str | None = dataclasses.field(default='none', metadata={
        'choices': ['none', 'very_strong', 'strong', 'moderate', 'weak', 'very_weak'],
        'show_default': True,
    })
    rates: str | None = dataclasses.field(default='U', metadata={
        'choices': ['G+I', 'G', 'I', 'U'],
        'show_default': True,

    })
    bootstraps: int | None = dataclasses.field(default=100, metadata={
        'help': 'Number of bootstrap replications',
        'show_default': True
    })
    ml_method: str | None = dataclasses.field(default='spr3', metadata={
        'choices': ['nni', 'spr3', 'spr5'],
        'show_default': True,
    })
    site_cov_cutoff: int | None = dataclasses.field(default=10, metadata={'show_default': True})
    threads: int | None = dataclasses.field(default=4, metadata={'show_default': True})
    include_ref: bool = dataclasses.field(default=False, metadata={'help': 'If set, reference is included in phylogeny'})
    include_filt_pos: bool = dataclasses.field(default=False)


class MainMega(BaseScript[Input, Output, Options]):
    """
    This class contains the main script for the MEGA model selection and tree building.
    """

    SNP_MATRIX_FILENAME = 'snp_matrix.fasta'

    def __init__(self, script_in: Input, script_out: Output, script_opts: Options) -> None:
        """
        Initializes the main script.
        :param script_in: Script input
        :param script_out: Script output
        ;param script_opts: Script options
        :return: None
        """
        tool_version = MLTreeConstruction().version
        super().__init__(
            name='MEGA',
            version=f'{tool_version}+CAMEL_{__VERSION__}',
            script_in=script_in,
            script_out=script_out,
            script_opts=script_opts
        )

    def _execute(self) -> None:
        """
        Runs the tool.
        :return: None
        """
        # Prepare SNP matrix
        if self._script_in.fasta is not None:
            fasta_path = self._script_opts.working_dir / MainMega.SNP_MATRIX_FILENAME
            fasta_path.symlink_to(self._script_in.fasta)
        else:
            vcf_files = [Path(v) for v, _ in self._script_in.vcf]
            sample_names = [n for _, n in self._script_in.vcf]
            fasta_path = self.__build_snp_matrix(
                vcf_files, sample_names, self._script_opts.include_ref, self._script_opts.include_filt_pos)
            if self._script_out.output_snp_matrix is not None:
                shutil.copyfile(fasta_path, self._script_out.output_snp_matrix)

        # Perform required actions
        if self._script_opts.action == 'model':
            self.__run_model_selection(fasta_path)
        elif self._script_opts.action == 'tree':
            self.__run_tree_construction(fasta_path, self._script_opts.model, self._script_opts.rates)
        else:
            model, rates = self.__run_model_selection(fasta_path)
            self.__run_tree_construction(fasta_path, model, rates)

    def __build_snp_matrix(self, vcf_files: list[Path], sample_names: list[str], include_ref: bool = False,
                           include_filtered_pos: bool = False) -> Path:
        """
        Builds a SNP matrix by combining all SNP positions across multiple VCF files.
        :param vcf_files: VCF files
        :param include_ref: If true, the reference is included in the phylogeny
        :param include_filtered_pos: If true, filtered positions are retained in the SNP matrix (as Ns)
        :return: SNP matrix path
        """
        snp_matrix_constructor = SnpMatrixConstructor()
        snp_matrix_constructor.update_parameters(
            include_ref=include_ref, include_filtered=None if include_filtered_pos else False)
        snp_matrix_constructor.add_input_files({
            'VCF': [ToolIOFile(v) for v in vcf_files],
            'SAMPLE_NAME': [ToolIOValue(n) for n in sample_names]})
        snp_matrix_constructor.run(self._script_opts.working_dir)
        return snp_matrix_constructor.tool_outputs['FASTA'][0].path

    def __run_model_selection(self, fasta_path: Path) -> tuple[str, str]:
        """
        Runs the model selection.
        :param fasta_path: FASTA input file
        :return: Model, rates parameter
        """
        model_selection = ModelSelection()
        MEGAUtils.update_model_selection_parameters(
            model_selection, self._script_opts.missing_data, self._script_opts.branch_swap, self._script_opts.site_cov_cutoff,
            self._script_opts.threads)
        model_selection.add_input_files({'FASTA': [ToolIOFile(fasta_path)]})
        dir_model_selection = self._script_opts.working_dir / 'model_selection'
        if not dir_model_selection.is_dir():
            dir_model_selection.mkdir()
        model_selection.run(dir_model_selection)
        shutil.copyfile(model_selection.tool_outputs['CSV'][0].path, self._script_out.output_model)
        return model_selection.informs['model'], model_selection.informs['rates_among_sites']

    def __run_tree_construction(self, fasta_path: Path, model: str, rates: str) -> None:
        """
        Runs the tree construction.
        :param fasta_path: FASTA input file
        :param model: Nucleotide substitution model
        :param rates: Rates among sites
        :return: None
        """
        tree_building = MLTreeConstruction()
        tree_building.add_input_files({'FASTA': [ToolIOFile(fasta_path)]})
        MEGAUtils.update_tree_building_parameters(
            tree_building, model, rates, self._script_opts.bootstraps, self._script_opts.missing_data, self._script_opts.site_cov_cutoff,
            self._script_opts.ml_method, self._script_opts.branch_swap, self._script_opts.threads)
        dir_tree_building = self._script_opts.working_dir / 'tree_building'
        if not dir_tree_building.is_dir():
            dir_tree_building.mkdir()
        tree_building.run(dir_tree_building)
        shutil.copyfile(tree_building.tool_outputs['NWK'][0].path, self._script_out.output_tree)

@click.command(name='mega', short_help='Model selection and/or tree construction using MEGA')
@click.option('--vcf', type=(click.Path(path_type=Path), str), multiple=True)
@cliutils.add_click_options_from_dataclass(Input, skip=['vcf'])
@cliutils.add_click_options_from_dataclass(Output)
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs: Any) -> None:
    """
    Runs the main script.
    :param kwargs: Command line arguments
    :return: None
    """
    main_script = MainMega(
        script_in = Input(**cliutils.from_kwargs(Input, kwargs)),
        script_out = Output(**cliutils.from_kwargs(Output, kwargs)),
        script_opts = Options(**cliutils.from_kwargs(Options, kwargs))
    )
    main_script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
