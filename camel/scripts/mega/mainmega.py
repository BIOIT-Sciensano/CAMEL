#!/usr/bin/env python
import argparse
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.scriptutils.basescript import BaseScript
from camel.app.loggers import initialize_logging
from camel.app.tools.mega.mltreeconstruction import MLTreeConstruction
from camel.app.tools.mega.modelselection import ModelSelection
from camel.app.tools.snpmatrix.snpmatrixconstructor import SnpMatrixConstructor
from camel.app.toolkits.phylogeny.megautils import MEGAUtils


class MainMega(BaseScript):
    """
    This class contains the main script for the MEGA model selection and tree building.
    """

    SNP_MATRIX_FILENAME = 'snp_matrix.fasta'

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        """
        super().__init__(name='MEGA', version='1.0', snakefile=None)
        self._args = MainMega._parse_arguments(args)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--fasta', type=Path, help="Input SNP matrix FASTA file")
        argument_parser.add_argument('--vcf', action='append', nargs=2)
        argument_parser.add_argument('--output-tree', type=Path, help="Output Newick tree file")
        argument_parser.add_argument('--output-model', type=Path, help="Output tabular model selection file")
        argument_parser.add_argument('--output-snp-matrix', type=Path,
                                     help='If set, the SNP matrix is exported to this FASTA file')
        argument_parser.add_argument('--action', choices=['both', 'model', 'tree'], required=True)
        argument_parser.add_argument('--missing-data', choices=[
            'use_all_sites', 'complete_deletion', 'partial_deletion'], default='use_all_sites')
        argument_parser.add_argument('--branch-swap', choices=[
            'none', 'very_strong', 'strong', 'moderate', 'weak', 'very_weak'], default='none')
        argument_parser.add_argument('--bootstraps', help="Number of bootstrap replications.",
                                     type=int, default=100)
        argument_parser.add_argument('--ml-method', choices=['nni', 'spr3', 'spr5'], default='spr3')
        argument_parser.add_argument('--site-cov-cutoff', choices=range(0, 101), type=int, default=50)
        argument_parser.add_argument('--model', choices=['JC', 'K2', 'T92', 'HKY', 'TN93', 'GTR'], default='JC')
        argument_parser.add_argument('--rates', choices=['G+I', 'G', 'I', 'U'], default='U')
        argument_parser.add_argument('--working-dir', type=Path, default=Path.cwd())
        argument_parser.add_argument(
            '--include-ref', action='store_true', help='If set, reference is included in phylogeny')
        argument_parser.add_argument(
            '--include-filt-pos', action='store_true',
            help='If True, filtered positions are retained in the SNP matrix as Ns')
        argument_parser.add_argument('--threads', type=int, default=3)
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the tool.
        :return: None
        """
        # Prepare SNP matrix
        if self._args.fasta is not None:
            fasta_path = self._args.working_dir / MainMega.SNP_MATRIX_FILENAME
            fasta_path.symlink_to(self._args.fasta)
        else:
            vcf_files = [Path(v) for v, _ in self._args.vcf]
            sample_names = [n for _, n in self._args.vcf]
            fasta_path = self.__build_snp_matrix(
                vcf_files, sample_names, self._args.include_ref, self._args.include_filt_pos)
            if self._args.output_snp_matrix is not None:
                shutil.copyfile(fasta_path, self._args.output_snp_matrix)

        # Perform required actions
        if self._args.action == 'model':
            self.__run_model_selection(fasta_path)
        elif self._args.action == 'tree':
            self.__run_tree_construction(fasta_path, self._args.model, self._args.rates)
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
        snp_matrix_constructor.run(self._args.working_dir)
        return snp_matrix_constructor.tool_outputs['FASTA'][0].path

    def __run_model_selection(self, fasta_path: Path) -> tuple[str, str]:
        """
        Runs the model selection.
        :param fasta_path: FASTA input file
        :return: Model, rates parameter
        """
        model_selection = ModelSelection()
        MEGAUtils.update_model_selection_parameters(
            model_selection, self._args.missing_data, self._args.branch_swap, self._args.site_cov_cutoff,
            self._args.threads)
        model_selection.add_input_files({'FASTA': [ToolIOFile(fasta_path)]})
        dir_model_selection = self._args.working_dir / 'model_selection'
        if not dir_model_selection.is_dir():
            dir_model_selection.mkdir()
        model_selection.run(dir_model_selection)
        shutil.copyfile(model_selection.tool_outputs['CSV'][0].path, self._args.output_model)
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
            tree_building, model, rates, self._args.bootstraps, self._args.missing_data, self._args.site_cov_cutoff,
            self._args.ml_method, self._args.branch_swap, self._args.threads)
        dir_tree_building = self._args.working_dir / 'tree_building'
        if not dir_tree_building.is_dir():
            dir_tree_building.mkdir()
        tree_building.run(dir_tree_building)
        shutil.copyfile(tree_building.tool_outputs['NWK'][0].path, self._args.output_tree)


if __name__ == '__main__':
    initialize_logging()
    main = MainMega()
    main.run()
