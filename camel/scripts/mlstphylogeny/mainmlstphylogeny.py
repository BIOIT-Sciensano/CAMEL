#!/usr/bin/env python
import argparse
import itertools
import shutil
import tempfile
from collections.abc import Sequence
from importlib.resources import files
from pathlib import Path
from typing import Optional

import pandas as pd

from camel.app.components import mainscriptutils
from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.phylogeny import mlstphyloutils
from camel.app.components.phylogeny.newickutils import NewickUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger, initialize_logging
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.figtree.figtree import FigTree
from camel.app.tools.grapetree.grapetree import GrapeTree


class MainMLSTPhylogeny:
    """
    Creates phylogenies based on the sequence typing output.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main scripts.
        :param args: Command line arguments
        """
        self._args = MainMLSTPhylogeny._parse_arguments(args)

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        report = self.__initialize_report()

        # Parse and filter the allele matrix
        allele_data = self.__parse_input_files()
        allele_data_filtered, cutoff_loci, cutoff_datasets = self.__filter_allele_matrix(allele_data)
        path_allele_matrix = self.__save_allele_matrix(allele_data_filtered)
        self.__report_add_section_filtering(
            report, path_allele_matrix, allele_data, allele_data_filtered, cutoff_loci, cutoff_datasets)

        # Calculate the distance matrix
        data_dist_matrix = self.__calculate_distance_matrix(allele_data_filtered)
        self.__report_add_section_distance_matrix(report, data_dist_matrix)

        # Check the distance matrix
        if all(max(values) == 0 for _, values in data_dist_matrix.iterrows()):
            self.__report_add_section_phylogeny_empty(report)
        else:
            # Create and visualize the phylogeny
            grapetree = self.run_grapetree(path_allele_matrix)

            # Store the tree image in a temporary file before adding it to the report
            with tempfile.NamedTemporaryFile(dir=str(self._args.dir_working), prefix='figtree_', suffix='.png') as \
                    nwk_out:
                figtree = NewickUtils.create_image_figtree(
                    grapetree.tool_outputs['NWK'][0].path,
                    Path(str(files('camel').joinpath('resources/figtree/template_cgmlst_tree.txt'))),
                    Path(nwk_out.name), self._args.image_width,
                    NewickUtils.calculate_tree_image_height(self._args.image_min_height, len(allele_data_filtered))
                )
                self.__report_add_section_phylogeny(report, grapetree, figtree)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: (optional) arguments
        :return: Arguments
        """
        ap = argparse.ArgumentParser()

        # Input files
        grp_input = ap.add_mutually_exclusive_group(required=True)
        grp_input.add_argument('--input-html', type=Path, action='append')
        grp_input.add_argument('--input-tsv', nargs=2, action='append')
        ap.add_argument('--html-key', help='Key for the scheme to use for the HTML input', default='cgmlst')
        ap.add_argument('--detection-method', type=str, choices=['blast', 'kma', 'rapid'], default='blast')
        ap.add_argument(
            '--tree-method', type=str, choices=['MSTreeV2', 'MSTree', 'NJ', 'RapidNJ', 'ninja', 'distance'],
            help='Tree building method for GrapeTree', default='MSTreeV2')

        # Allele matrix filtering
        ap.add_argument(
            '--min-perc-loci', type=int, default=90,
            help='Minimum percentage of loci that should be present in a dataset')
        ap.add_argument(
            '--min-perc-samples', type=int, default=90,
            help='Minimum percentage of datasets where loci should be present')
        ap.add_argument(
            '--keep-all-loci', action='store_true',
            help='If true, all loci are retained in the allele matrix output (including loci absent in all datasets')
        ap.add_argument('--no-temp-allele-ids', action='store_true',
                        help='If enabled, the temporary (hashed) allele ids are not used')

        # Output files
        ap.add_argument('--dir-working', type=Path, help='Working directory', default=Path.cwd())
        ap.add_argument('--output-allele-matrix', type=Path, help='Output path for allele matrix')
        ap.add_argument('--output-dist-matrix', type=Path, help='Output path for the AD matrix')
        ap.add_argument('--output-tree', type=Path, help='Output path for the Newick tree')
        ap.add_argument('--output-html', type=Path, help='Path to the HTML output report', required=True)
        ap.add_argument('--output-dir', type=Path, help='Path to the output directory', required=True)

        # Visualization
        ap.add_argument(
            '--image-min-height', type=int, default=256,
            help="Minimum height for image, is increased based on the number of isolates that are plotted")
        ap.add_argument('--image-width', type=int, default=640)
        return ap.parse_args(args)

    @staticmethod
    def calculate_distance(alleles_a: pd.Series, alleles_b: pd.Series) -> int:
        """
        Calculates the number of different alleles.
        :param alleles_a: Alleles a
        :param alleles_b: Alleles b
        :return: Number of allelic differences
        """
        total_dist = 0
        for allele_a, allele_b in zip(alleles_a, alleles_b):
            if allele_a == '-' and allele_b == '-':
                continue
            else:
                total_dist += 0 if allele_a == allele_b else 1
        return total_dist

    def __initialize_report(self) -> HtmlReport:
        """
        Initializes the HTML report
        :return: HTML report
        """
        # Header
        report = mainscriptutils.init_report(
            self._args.output_html, self._args.output_dir, 'MLST phylogeny report', 'MLST phylogeny')
        if self._args.input_tsv is not None:
            input_file_str = ', '.join([Path(name).name for _, name in self._args.input_tsv])
        else:
            input_file_str = f'Sequence typing output ({len(self._args.input_html)} datasets)'

        # Analysis info section
        report.add_html_object(mainscriptutils.generate_analysis_info_section(
            self._args, additional_info=[
                ['Allele detection method:', self._args.detection_method],
                ['Tree building method:', self._args.tree_method]
            ],
            input_file_str=input_file_str))
        report.save()
        return report

    def __parse_input_files(self) -> pd.DataFrame:
        """
        Parses the input files for the tree construction.
        :return: Detected alleles by dataset name
        """
        if self._args.input_html:
            allele_data = mlstphyloutils.parse_html_typing_list(
                self._args.input_html, self._args.html_key, self._args.detection_method,
                not self._args.no_temp_allele_ids)
        elif self._args.input_tsv:
            allele_data = mlstphyloutils.parse_tsv_typing_list(
                [(Path(path), name) for path, name in self._args.input_tsv], self._args.detection_method,
                not self._args.no_temp_allele_ids)
        else:
            raise ValueError("No input files specified")
        if len(allele_data) < 3:
            raise ValueError("At least 3 datasets are required")
        logger.info(f"Alleles parsed for {len(allele_data)} input files ({len(allele_data.columns)} loci)")
        return allele_data

    def __filter_allele_matrix(self, allele_data: pd.DataFrame) -> tuple[pd.DataFrame, int, int]:
        """
        Filters the allele matrix by removing:
        - Datasets with less than x% of loci detected
        - Loci present in less than x% of datasets
        :param allele_data: Allele data
        :return: Filtered allele data, loci cutoff, datasets cutoff
        """
        # Filter allele matrix (nb. of loci detected per dataset)
        nb_loci_detected = allele_data.apply(lambda x: len(x) - list(x).count('-'), axis=1)
        cutoff_loci = int(self._args.min_perc_loci * len(allele_data.columns) / 100)
        logger.info(f"Removing datasets with < {cutoff_loci} ({self._args.min_perc_loci}%) loci detected")
        allele_data_filt = allele_data[nb_loci_detected > cutoff_loci]
        logger.info(f"{len(allele_data_filt)} datasets passed filtering")

        # Filter allele matrix (loci detected in nb. of datasets)
        locus_present_in_datasets = allele_data_filt.apply(lambda x: len(x) - list(x).count('-'))
        cutoff_datasets = int(self._args.min_perc_samples * len(allele_data_filt) / 100)
        if not self._args.keep_all_loci:
            logger.info(f"Removing loci detected < {cutoff_datasets} ({self._args.min_perc_samples}%) datasets")
            allele_data_filt = allele_data_filt.iloc[:, list(locus_present_in_datasets > cutoff_datasets)]
        logger.info(f"{len(allele_data_filt.columns)} loci passed filtering")
        return allele_data_filt, cutoff_loci, cutoff_datasets

    def __save_allele_matrix(self, allele_data_filt: pd.DataFrame) -> Path:
        """
        Saves the allele matrix to a file.
        :param allele_data_filt: Filtered allele data
        """
        path_allele_matrix = self._args.dir_working / 'allele_matrix-filtered.tsv'
        allele_data_filt.to_csv(path_allele_matrix, index_label='ID', sep='\t')

        # Save the allele matrix to a separate file is specified
        if self._args.output_allele_matrix:
            shutil.copyfile(path_allele_matrix, self._args.output_allele_matrix)
            logger.info(f"Filtered allele matrix saved to: {self._args.output_allele_matrix}")
        return path_allele_matrix

    def __report_add_section_filtering(
            self, report: HtmlReport, path_allele_matrix: Path, allele_data: pd.DataFrame,
            allele_data_filt: pd.DataFrame, cutoff_loci: int, cutoff_datasets: int) -> Path:
        """
        Adds the filtering section to the output report.
        :param report: HTML report
        :param path_allele_matrix: Path to the allele matrix
        :param allele_data: Original allele data
        :param allele_data_filt: Filtered allele data
        :param cutoff_loci: Nb. of loci cutoff value
        :param cutoff_datasets: Nb. of datasets cutoff value
        :return: Path to the allele matrix
        """
        section_filtering = HtmlReportSection('Allele matrix filtering')
        section_filtering.add_table([
            ['Loci required (%):', f'{self._args.min_perc_loci}%'],
            ['Loci required:', str(cutoff_loci)],
            ['Present in datasets (%):', f'{self._args.min_perc_samples}%'],
            ['Present in datasets', str(cutoff_datasets)]
        ], None, [('class', 'information')])
        section_filtering.add_table([
            ['Before filtering', len(allele_data), len(allele_data.columns)],
            ['After filtering', len(allele_data_filt), len(allele_data_filt.columns)]
        ], ['Step', 'Datasets', 'Loci'], [('class', 'data')])
        relative_path = Path('allele_matrix-filtered.tsv')
        section_filtering.add_link_to_file('Download filtered allele matrix (TSV)', relative_path)
        section_filtering.add_file(path_allele_matrix, relative_path)
        section_filtering.add_warning_message('Locus filtering is performed after low quality datasets are removed.')
        report.add_html_object(section_filtering)
        section_filtering.copy_files(report.output_dir)
        report.save()
        return path_allele_matrix

    def __calculate_distance_matrix(self, allele_data_filtered: pd.DataFrame) -> pd.DataFrame:
        """
        Calculates the pairwise distance matrix
        :param allele_data_filtered: Filtered allele matrix
        :return: Distance matrix
        """
        # Calculate pair-wise distances
        distance_by_dataset_pair = {}
        for dataset_a, dataset_b in itertools.combinations(allele_data_filtered.index, r=2):
            key = tuple(sorted([dataset_a, dataset_b]))
            dist = MainMLSTPhylogeny.calculate_distance(
                allele_data_filtered.loc[dataset_a], allele_data_filtered.loc[dataset_b])
            distance_by_dataset_pair[key] = dist

        # Create data frame with pairwise distances
        records_out = []
        for dataset_a in allele_data_filtered.index:
            records_out.append({
                dataset_b: distance_by_dataset_pair.get(tuple(sorted([dataset_a, dataset_b])), 0) for
                dataset_b in allele_data_filtered.index
            })
        return pd.DataFrame(records_out, index=allele_data_filtered.index)

    def __report_add_section_distance_matrix(self, report: HtmlReport, data_dist_matrix: pd.DataFrame) -> None:
        """
        Adds the distance matrix section to the HTML report.
        :param report: HTML report
        :param data_dist_matrix: Distance matrix data
        :return: None
        """
        # Save distance matrix in output folder
        path_tsv_dist = self._args.dir_working / 'dist_matrix.tsv'
        data_dist_matrix.to_csv(path_tsv_dist, sep='\t', index_label='ID')
        # Copy the distance file (if specified)
        if self._args.output_dist_matrix is not None:
            shutil.copyfile(path_tsv_dist, self._args.output_dist_matrix)
            logger.info(f"Distance matrix saved to: {self._args.output_dist_matrix}")

        # Section with pair-wise distances
        section_distances = HtmlReportSection('Allele distances (AD)')
        header = [HtmlElement('th', '')] + [
            HtmlElement('th', col, [('style', 'writing-mode: vertical-lr;')]) for col in data_dist_matrix.columns]
        # noinspection PyTypeChecker
        section_distances.add_table(
            [header] + [[str(dataset)] + [str(x) for x in data] for dataset, data in data_dist_matrix.iterrows()],
            None, [('class', 'data')])
        relative_path = Path('allele_dist_matrix.tsv')
        section_distances.add_file(path_tsv_dist, relative_path)
        section_distances.add_link_to_file('Download AD matrix (TSV)', relative_path)
        section_distances.add_warning_message('Allele distances are calculated based on the filtered allele matrix.')
        report.add_html_object(section_distances)
        section_distances.copy_files(report.output_dir)
        report.save()

    def run_grapetree(self, allele_matrix: Path) -> GrapeTree:
        """
        Runs GrapeTree to create the phylogeny.
        :param allele_matrix: Input allele matrix
        :return: GrapeTree tool instance
        """
        dir_grapetree = self._args.dir_working / 'grapetree'
        dir_grapetree.mkdir(exist_ok=True, parents=True)
        grapetree = GrapeTree()
        output_path = dir_grapetree / 'tree.nwk'
        grapetree.update_parameters(output_path=str(output_path))
        grapetree.add_input_files({'TSV': [ToolIOFile(allele_matrix)]})
        grapetree.run(dir_grapetree)
        logger.info(f"Phylogeny created: {output_path}")

        # Copy the Newick file (if specified)
        path_newick = grapetree.tool_outputs['NWK'][0].path
        if self._args.output_tree is not None:
            shutil.copyfile(path_newick, self._args.output_tree)
            logger.info(f"Newick saved to: {self._args.output_tree}")
        return grapetree

    def __report_add_section_phylogeny(self, report: HtmlReport, grapetree: GrapeTree, figtree: FigTree) -> None:
        """
        Adds the phylogeny section to the HTML report.
        :param report: HTML report
        :param grapetree: GrapeTree tool object
        :param figtree: FigTree tool object
        :return: None
        """
        # Add output tree to the output
        section_phylo = HtmlReportSection('Phylogeny')
        relative_path_img = Path('phylogeny.png')
        section_phylo.add_file(figtree.tool_outputs['PNG'][0].path, relative_path_img)
        section_phylo.add_html_object(HtmlElement('img', attributes=[('src', str(relative_path_img)), ('border', '1')]))
        section_phylo.add_line_break()
        relative_path = Path('phylogeny.nwk')
        section_phylo.add_link_to_file('Download tree (NWK)', relative_path)
        section_phylo.add_file(grapetree.tool_outputs['NWK'][0].path, relative_path)
        section_phylo.copy_files(report.output_dir)
        report.add_html_object(section_phylo)

        # Add commands and citations
        report.add_html_object(SnakePipelineUtils.create_commands_section(
            [grapetree.informs, figtree.informs], self._args.dir_working))
        report.add_html_object(SnakePipelineUtils.create_citations_section(['Zhou_2018-grapetree']))
        report.save()

    def __report_add_section_phylogeny_empty(self, report: HtmlReport) -> None:
        """
        Adds the phylogeny section to the HTML report.
        :param report: HTML report
        """
        section = HtmlReportSection('Phylogeny')
        section.add_error_message('All profiles are identical, cannot generate phylogeny.')
        report.add_html_object(section)
        report.save()


def main(args: Sequence[str] | None = None) -> None:
    """
    Entry point for the common interface.
    :param args: Command line arguments
    :return: None
    """
    script = MainMLSTPhylogeny(args)
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
