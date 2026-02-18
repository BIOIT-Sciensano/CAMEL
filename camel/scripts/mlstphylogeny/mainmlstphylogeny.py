#!/usr/bin/env python
import dataclasses
import shutil
import tempfile
from importlib.resources import files
from pathlib import Path

import click
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

from camel.app.cli import cliutils
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.reports import reportutils
from camel.app.core.reports.htmlelement import HtmlElement
from camel.app.core.reports.htmlreport import HtmlReport
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.loggers import logger, initialize_logging
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.model import BaseOptions, BaseOutput, BaseInput
from camel.app.toolkits.phylogeny import mlstphyloutils
from camel.app.toolkits.phylogeny.newickutils import NewickUtils
from camel.app.tools.figtree.figtree import FigTree
from camel.app.tools.grapetree.grapetree import GrapeTree


@dataclasses.dataclass(frozen=True)
class Input(BaseInput):
    """
    Defines the script input.
    """

    input_tsv: list[tuple[Path, str]] | None
    input_html: list[Path] | None
    html_key: str | None = None
    detection_method: str = dataclasses.field(default='blast', metadata={"choices": ["blast", "mist", "kma"]})

@dataclasses.dataclass(frozen=True)
class Output(BaseOutput):
    """
    Defines the script output.
    """

    output_html: Path = None
    output_dir: Path | None = None
    output_allele_matrix: Path | None = None
    output_dist_matrix: Path | None = None
    output_tree: Path | None = None

@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Defines the script options.
    """

    dir_working: Path = dataclasses.field(default=Path.cwd(), metadata={"help": "Working directory"})
    tree_method: str = dataclasses.field(default="MSTreeV2", metadata={"help": "Tree building method"})
    min_perc_loci: int = dataclasses.field(default=90, metadata={"help": "Minimum percentage of loci per dataset"})
    min_perc_samples: int = dataclasses.field(default=90, metadata={"help": "Minimum percentage of datasets where loci should be present"},)
    keep_all_loci: bool = dataclasses.field(default=False, metadata={"help": "Retain all loci in allele matrix"})
    no_temp_allele_ids: bool = dataclasses.field(default=False, metadata={"help": "Do not use temporary hashed allele IDs"})
    image_width: int = dataclasses.field(default=640, metadata={"help": "Width of tree image"})
    image_min_height: int = dataclasses.field(default=256, metadata={"help": "Minimum height of tree image"})


class MainMLSTPhylogeny(BaseScript[Input, Output, Options]):
    """
    Creates phylogenies based on the sequence typing output.
    """

    def __init__(self, in_: Input, out: Output, opts: Options) -> None:
        """
        Initializes the main scripts.
        :param in_: Script input
        :param out: Script output
        :param opts: Script options
        :return: None
        """
        super().__init__(
            name='MLST phylogeny',
            version='1.0',
            script_in=in_,
            script_out=out,
            script_opts=opts
        )

    def _execute(self) -> None:
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
        data_dist_matrix = self.pdist_calc_distance_matrix(allele_data_filtered)
        self.__report_add_section_distance_matrix(report, data_dist_matrix)

        # Check the distance matrix
        if all(max(values) == 0 for _, values in data_dist_matrix.iterrows()):
            self.__report_add_section_phylogeny_empty(report)
        else:
            # Create and visualize the phylogeny
            grapetree = self.run_grapetree(path_allele_matrix)

            # Store the tree image in a temporary file before adding it to the report
            with tempfile.NamedTemporaryFile(dir=str(self._script_opts.dir_working), prefix='figtree_', suffix='.png') as \
                    nwk_out:
                figtree = NewickUtils.create_image_figtree(
                    grapetree.tool_outputs['NWK'][0].path,
                    Path(str(files('camel').joinpath('resources/tools/figtree/template_cgmlst_tree.txt'))),
                    Path(nwk_out.name), self._script_opts.image_width,
                    NewickUtils.calculate_tree_image_height(self._script_opts.image_min_height, len(allele_data_filtered))
                )
                self.__report_add_section_phylogeny(report, grapetree, figtree)

    @staticmethod
    def encode_alleles_global(df: pd.DataFrame) -> pd.DataFrame:
        """
        Converts allele ids from string to integers for faster comparison.
        :param df: Input dataframe
        :return: Encoded dataframe
        """
        df = df.replace('-', -1)
        flat = pd.factorize(df.values.ravel())[0]
        return pd.DataFrame(
            flat.reshape(df.shape), index=df.index, columns=df.columns, dtype=np.int32
        )

    @staticmethod
    def pdist_allele_distance(u: np.ndarray, v: np.ndarray) -> np.int32:
        """
        Calculates the allele distance between two arrays.
        :param u: Array A (allele calls for a dataset)
        :param v: Array B (allele calls for a dataset)
        :return: The total distance
        """
        both_missing = (u == -1) & (v == -1)
        diff = u != v
        # noinspection PyUnresolvedReferences
        diff[both_missing] = False
        return np.sum(diff)

    @staticmethod
    def pdist_calc_distance_matrix(df_alleles: pd.DataFrame) -> pd.DataFrame:
        """
        Uses pdist to calculate the pairwise distances.
        :param df_alleles: Input dataframe with alleles (encoded as integers)
        :return: Pairwise distance matrix
        """
        dists = pdist(df_alleles.to_numpy(), metric=MainMLSTPhylogeny.pdist_allele_distance)
        return pd.DataFrame(
            squareform(dists),
            index=df_alleles.index,
            columns=df_alleles.index,
            dtype=int
        )

    def __initialize_report(self) -> HtmlReport:
        """
        Initializes the HTML report.
        :return: HTML report
        """
        # Header
        report = reportutils.init_report(
            path_out=self._script_out.output_html,
            key='MLST phylogeny',
            title='MLST phylogeny report',
            dir_out=self._script_out.output_dir)
        if self._script_in.input_tsv is not None:
            input_file_str = f'Sequence typing output ({len(self._script_in.input_tsv)} datasets)'
        else:
            input_file_str = f'Sequence typing output ({len(self._script_in.input_html)} datasets)'

        # Analysis info section
        report.add_html_object(reportutils.create_overview_section(
            version=self._version,
            dataset_name=input_file_str,
            input_file_str=input_file_str,
            extra_data=[
                ['Allele detection method', self._script_in.detection_method],
                ['Tree building method', self._script_opts.tree_method]
            ]
        ))
        report.save()
        return report

    def __parse_input_files(self) -> pd.DataFrame:
        """
        Parses the input files for the tree construction.
        :return: Detected alleles by dataset name
        """
        if self._script_in.input_html:
            allele_data = mlstphyloutils.parse_html_typing_list(
                self._script_in.input_html, self._script_in.html_key, self._script_in.detection_method,
                not self._script_opts.no_temp_allele_ids)
        elif self._script_in.input_tsv:
            allele_data = mlstphyloutils.parse_tsv_typing_list(
                [(Path(path), name) for path, name in self._script_in.input_tsv], self._script_in.detection_method,
                not self._script_opts.no_temp_allele_ids)
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
        cutoff_loci = int(self._script_opts.min_perc_loci * len(allele_data.columns) / 100)
        logger.info(f"Removing datasets with < {cutoff_loci} ({self._script_opts.min_perc_loci}%) loci detected")
        allele_data_filt = allele_data[nb_loci_detected > cutoff_loci]
        logger.info(f"{len(allele_data_filt)} datasets passed filtering")

        # Filter allele matrix (loci detected in nb. of datasets)
        locus_present_in_datasets = allele_data_filt.apply(lambda x: len(x) - list(x).count('-'))
        cutoff_datasets = int(self._script_opts.min_perc_samples * len(allele_data_filt) / 100)
        if not self._script_opts.keep_all_loci:
            logger.info(f"Removing loci detected < {cutoff_datasets} ({self._script_opts.min_perc_samples}%) datasets")
            allele_data_filt = allele_data_filt.iloc[:, list(locus_present_in_datasets > cutoff_datasets)]
        logger.info(f"{len(allele_data_filt.columns)} loci passed filtering")
        return allele_data_filt, cutoff_loci, cutoff_datasets

    def __save_allele_matrix(self, allele_data_filt: pd.DataFrame) -> Path:
        """
        Saves the allele matrix to a file.
        :param allele_data_filt: Filtered allele data
        """
        path_allele_matrix = self._script_opts.dir_working / 'allele_matrix-filtered.tsv'
        allele_data_filt.to_csv(path_allele_matrix, index_label='ID', sep='\t')

        # Save the allele matrix to a separate file is specified
        if self._script_out.output_allele_matrix:
            shutil.copyfile(path_allele_matrix, self._script_out.output_allele_matrix)
            logger.info(f"Filtered allele matrix saved to: {self._script_out.output_allele_matrix}")
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
            ['Loci required (%):', f'{self._script_opts.min_perc_loci}%'],
            ['Loci required:', str(cutoff_loci)],
            ['Present in datasets (%):', f'{self._script_opts.min_perc_samples}%'],
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

    def __report_add_section_distance_matrix(self, report: HtmlReport, data_dist_matrix: pd.DataFrame) -> None:
        """
        Adds the distance matrix section to the HTML report.
        :param report: HTML report
        :param data_dist_matrix: Distance matrix data
        :return: None
        """
        # Save distance matrix in output folder
        path_tsv_dist = self._script_opts.dir_working / 'dist_matrix.tsv'
        data_dist_matrix.to_csv(path_tsv_dist, sep='\t', index_label='ID')
        # Copy the distance file (if specified)
        if self._script_out.output_dist_matrix is not None:
            shutil.copyfile(path_tsv_dist, self._script_out.output_dist_matrix)
            logger.info(f"Distance matrix saved to: {self._script_out.output_dist_matrix}")

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
        dir_grapetree = self._script_opts.dir_working / 'grapetree'
        dir_grapetree.mkdir(exist_ok=True, parents=True)
        grapetree = GrapeTree()
        output_path = dir_grapetree / 'tree.nwk'
        grapetree.update_parameters(output_path=str(output_path))
        grapetree.add_input_files({'TSV': [ToolIOFile(allele_matrix)]})
        grapetree.run(dir_grapetree)
        logger.info(f"Phylogeny created: {output_path}")

        # Copy the Newick file (if specified)
        path_newick = grapetree.tool_outputs['NWK'][0].path
        if self._script_out.output_tree is not None:
            shutil.copyfile(path_newick, self._script_out.output_tree)
            logger.info(f"Newick saved to: {self._script_out.output_tree}")
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
        report.add_html_object(reportutils.create_commands_section(
            [grapetree.informs, figtree.informs], self._script_opts.dir_working))
        report.add_html_object(reportutils.create_citations_section(['Zhou_2018-grapetree']))
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


@click.command(name='mlst_phylogeny', short_help='Phylogenetic investigation based on (cg)MLST outputs')
@cliutils.add_click_options_from_dataclass(Input, skip=['input_html', 'input_tsv'])
@click.option('--input-html', type=click.Path(exists=True, path_type=Path), multiple=True, help="Input HTML files")
@click.option('--input-tsv', nargs=2, type=(click.Path(exists=True, path_type=Path), str), multiple=True, help="Input TSV files (two per entry)")
@cliutils.add_click_options_from_dataclass(Output)
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Entry point for the common interface.
    :param kwargs: Command line arguments
    :return: None
    """
    script_in = Input(**cliutils.from_kwargs(Input, kwargs))
    script_out = Output(**cliutils.from_kwargs(Output, kwargs))
    script_opts = Options(**cliutils.from_kwargs(Options, kwargs))
    mlst_phylo = MainMLSTPhylogeny(in_=script_in, out=script_out, opts=script_opts)
    mlst_phylo.run()


if __name__ == '__main__':
    initialize_logging()
    main()
