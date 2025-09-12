import argparse
import datetime
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Optional, Any, Union

import pandas as pd
from Bio import SeqIO

from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.components.phylogeny.megautils import MEGAUtils
from camel.app.components.phylogeny.newickutils import NewickUtils
from camel.app.components.workflows.trimmingilluminawrapper import TrimmingIlluminaWrapper
from camel.app.config import config
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.snakemake import snakemakeutils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.mega.mltreeconstruction import MLTreeConstruction
from camel.app.tools.mega.modelselection import ModelSelection
from camel.app.tools.snpmatrix.snpmatrixconstructor import SnpMatrixConstructor
from camel.resources import CSS_STYLE
from camel.scripts.snpphylogeny import SNAKEFILE_TRIMMING_ALL
from camel.scripts.snpphylogeny.snakefile.trimming_all import TRIMMING_ALL


@dataclass(frozen=True)
class Sample:
    """
    This class is used to represent an input sample.
    """
    name_full: str
    name_valid: str
    reads_raw: list[ToolIOFile]
    reads_names: list[str]

    def __hash__(self) -> int:
        """
        Returns the hash value for this
        :return: Hash value
        """
        return hash(self.name_valid)

    def __eq__(self, other: 'Sample') -> bool:
        """
        Checks if two sample instances are equal.
        :param other: Other instance
        :return: True if equal, False otherwise
        """
        return self.name_valid == other.name_valid


@dataclass
class MappingInput:
    """
    This class contains the input for the read mapping step.
    """
    pe: list[ToolIOFile]
    se_fwd: Optional[ToolIOFile] = None
    se_rev: Optional[ToolIOFile] = None

    def as_dict(self) -> dict[str, Any]:
        """
        Return the mapping input as a dictionary
        :return: Mapping input in YAML format
        """
        mapping_data: dict[str, Union[list[str], str]] = {'PE': [str(x.path) for x in self.pe]}
        if self.se_fwd is not None:
            mapping_data['SE_FWD'] = str(self.se_fwd.path)
        if self.se_rev is not None:
            mapping_data['SE_REV'] = str(self.se_rev.path)
        return mapping_data


class SnpPhylogenyUtils:
    """
    Utility class for the SNP phylogeny pipelines.
    """

    @staticmethod
    def initialize_report(pipeline_name: str, args: argparse.Namespace) -> HtmlReport:
        """
        Initializes a HTML report.
        :param pipeline_name: Pipeline name (e.g. 'Samtools', 'PHEnix')
        :param args: Pipeline arguments
        :return: Initializes HTML report
        """
        output_dir = Path(args.output_dir)
        report = HtmlReport(Path(args.output_html), output_dir)
        if not output_dir.exists():
            output_dir.mkdir(parents=True)
        report.initialize(f'SNP Phylogeny - {pipeline_name}', CSS_STYLE)
        report.add_pipeline_header(f'SNP Phylogeny ({pipeline_name})')
        section = HtmlReportSection('Analysis info')
        section.add_table([
            ['Analysis date:', datetime.datetime.now().strftime(config.date_fmt)],
            ['Nb. of samples:', len(args.sample)],
            ['Reference:', args.reference_name if args.reference_name else Path(args.reference).name]],
            table_attributes=[('class', 'information')]
        )
        report.add_html_object(section)
        report.save()
        return report

    @staticmethod
    def add_common_arguments(argument_parser: argparse.ArgumentParser) -> None:
        """
        Adds the common arguments to the given argument parser.
        :param argument_parser: Argument parser
        :return: None
        """
        argument_parser.add_argument('--output-dir', required=True, type=Path)
        argument_parser.add_argument('--output-html', required=True, type=Path)
        argument_parser.add_argument('--working-dir', type=Path, default=Path.cwd())
        argument_parser.add_argument('--output-fasta', type=Path, help='SNP matrix output')
        argument_parser.add_argument('--threads', default=8, type=int)
        argument_parser.add_argument('--reference', required=True, type=Path)
        argument_parser.add_argument('--reference-name', type=str)
        argument_parser.add_argument('--sample', nargs=5, action='append', required=True)
        argument_parser.add_argument('--trim-reads', action='store_true')
        argument_parser.add_argument(
            '--adapter', help="(Illumina) Adapter that was used for sequencing, used for read-trimming",
            choices=['NexteraPE', 'TruSeq2', 'TruSeq3'], default='NexteraPE')
        argument_parser.add_argument('--missing-data', choices=[
            'complete_deletion', 'use_all_sites', 'partial_deletion'], default='use_all_sites')
        argument_parser.add_argument('--site-cov-cutoff', choices=range(0, 101), type=int)
        argument_parser.add_argument('--branch-swap', choices=[
            'none', 'weak', 'very_weak', 'moderate', 'strong', 'very_strong'], default='moderate')
        argument_parser.add_argument('--bootstraps', type=int, default=100)
        argument_parser.add_argument('--ml-method', choices=['nni', 'spr3', 'spr5'], default='spr3')

    @staticmethod
    def extract_samples(args: argparse.Namespace) -> list[Sample]:
        """
        Returns the sample names (sorted alphabetically).
        :param args: Command line arguments
        :return: Samples
        """
        samples = []
        for sample_name, fwd_name, fwd_read, rev_name, rev_read in args.sample:
            if (sample_name is None) or (len(FileSystemHelper.make_valid(sample_name)) == 0):
                sample_name = FastqUtils.get_sample_name(fwd_name)
            samples.append(Sample(
                sample_name, FileSystemHelper.make_valid(sample_name),
                [ToolIOFile(Path(fwd_read)), ToolIOFile(Path(rev_read))],
                [fwd_name, rev_name]
            ))

        # Check for duplicate names
        if len(set(s.name_valid for s in samples)) != len(args.sample):
            sample_names = [s.name_valid for s in samples]
            duplicate_sample_names = [x for x in set(sample_names) if sample_names.count(x) > 1]
            raise ValueError("Duplicate sample names are not allowed. Conflicting sample(s): {}".format(
                ', '.join([f"'{n}' ({sample_names.count(n)} times)" for n in duplicate_sample_names])))

        # Check for empty sample names
        empty_samples = [s for s in samples if len(s.name_valid) == 0]
        if len(empty_samples) > 0:
            raise ValueError("Empty sample name for sample(s): {}".format(', '.join([
                str(s.reads_names) for s in empty_samples])))

        # Sort samples by name
        return sorted(samples, key=lambda s: s.name_valid)

    @staticmethod
    def symlink_input_files(samples: list[Sample], dir_working: Path) -> dict[Sample, list[ToolIOFile]]:
        """
        This function creates symlinks for the raw read files that belong to the given samples.
        :param samples: Samples
        :param dir_working: Working directory
        :return: Symlink locations by sample
        """
        logger.info(f"Creating symlinks for input files for {len(samples)} samples")
        symlink_dir = dir_working / 'input'
        if not symlink_dir.exists():
            symlink_dir.mkdir(parents=True)
        fq_by_sample = {s: [] for s in samples}
        for sample in samples:
            for nb, fq_file in enumerate(sample.reads_raw, 1):
                is_gzipped = FileSystemHelper.is_gzipped(fq_file.path)
                path_link = symlink_dir / f"{sample.name_valid}_{nb}.fastq{'.gz' if is_gzipped else ''}"
                if path_link.exists():
                    path_link.unlink()
                path_link.symlink_to(fq_file.path)
                fq_by_sample[sample].append(ToolIOFile(path_link))
        return fq_by_sample

    @staticmethod
    def trim_all_reads(fq_by_sample: dict[Sample, list[ToolIOFile]], working_dir: Path, adapter: str,
                       threads: int = 8) -> dict[Sample, TrimmingIlluminaWrapper.ReadTrimmingOutput]:
        """
        Trims all the reads in parallel using Snakemake.
        :param fq_by_sample: FASTQ files by sample
        :param working_dir: Working directory
        :param adapter: Adapters to trim
        :param threads: Number of threads
        :return: None
        """
        config_data = {
            'working_dir': str(working_dir),
            'samples': {s.name_valid: [str(f.path) for f in fq] for s, fq in fq_by_sample.items()},
            'read_trimming': {'adapter': adapter}
        }
        config_file = SnakePipelineUtils.generate_config_file(config_data, working_dir)
        output_file = TRIMMING_ALL
        SnakePipelineUtils.run_snakemake(
            SNAKEFILE_TRIMMING_ALL, config_file, [output_file], working_dir, threads)
        trimming_out_by_sample = snakemakeutils.load_object(working_dir / output_file)
        return {s: trimming_out_by_sample[s.name_valid] for s in fq_by_sample}

    @staticmethod
    def add_trimming_section_empty(report: HtmlReport) -> None:
        """
        Adds an empty trimming section to the report.
        :param report: Report
        :return: None
        """
        section = HtmlReportSection('Read trimming')
        section.add_paragraph('Read trimming disabled.')
        report.add_html_object(section)
        report.save()

    @staticmethod
    def add_trimming_section(report: HtmlReport, trimming_output_by_sample: dict[
            Sample, TrimmingIlluminaWrapper.ReadTrimmingOutput]) -> None:
        """
        Adds the trimming section to the report.
        :param report: HTML report
        :param trimming_output_by_sample: Trimming output dictionary
        :return: None
        """
        section = HtmlReportSection('Read trimming')

        table_data = []
        for sample, trimming_output in trimming_output_by_sample.items():
            # Add stats
            stats = trimming_output.informs_trimming
            row = [sample.name_full, stats['paired_reads_in'], stats['paired_reads_out'], stats['forward_only_reads'],
                   stats['reverse_only_reads'], stats['reads_drop']]

            # Add FastQC reports
            for i, f in enumerate(trimming_output.fastq_reports_post, start=1):
                relative_path = Path('fastqc_report') / f'{sample.name_valid}_{i}.html'
                section.add_file(f.path, relative_path)
                row.append(HtmlTableCell('view', link=str(relative_path)))
            table_data.append(row)

        # Add table
        header = ['Sample', 'Total read pairs', 'Both surviving', 'Forward only surviving', 'Reverse only surviving',
                  'Dropped', 'FastQC report (fwd.)', 'FastQC report (rev.)']
        section.add_table(table_data, header, [('class', 'data')])

        # Save the report
        report.add_html_object(section)
        section.copy_files(report.output_dir)
        report.save()

    @staticmethod
    def construct_snp_matrix(sample_names: list[str], vcf_files: list[ToolIOFile], working_dir: Path,
                             include_ref: bool = False, include_filtered: bool = True) -> ToolIOFile:
        """
        Constructs a SNP matrix based on the given VCF files.
        :param sample_names: Sample names
        :param vcf_files: VCF files
        :param working_dir: Working directory
        :param include_ref: If True, the reference is included in the SNP matrix
        :param include_filtered: If True, filtered variants are included in the matrix
        :return: SNP matrix ToolIOFile
        """
        if not working_dir.exists():
            working_dir.mkdir(parents=True)
        snp_matrix_constructor = SnpMatrixConstructor()

        # Parameters and input
        if include_ref is True:
            snp_matrix_constructor.update_parameters(include_ref=None)
        if include_filtered is True:
            snp_matrix_constructor.update_parameters(include_filtered=None)
        snp_matrix_constructor.add_input_files({
            'VCF': vcf_files,
            'SAMPLE_NAME': [ToolIOValue(s) for s in sample_names],
        })

        # Run tool
        snp_matrix_constructor.run(working_dir)
        return snp_matrix_constructor.tool_outputs['FASTA'][0]

    @staticmethod
    def add_output_files_section(report: HtmlReport, column_names: list[str], output_files: dict[
            Sample, list[Optional[Path]]], snp_matrix: Path) -> None:
        """
        Adds the section with the output files.
        :param report: Report
        :param column_names: Output file column names
        :param output_files: Output files
        :param snp_matrix: Path to the SNP matrix
        :return: None
        """
        section = HtmlReportSection('Output files')

        # Add SNP matrix
        relative_path = Path('snp_matrix.fasta')
        section.add_file(snp_matrix, relative_path)
        table_data_snp_matrix = [
            ['SNP matrix:', HtmlTableCell('Download (FASTA)', link=str(relative_path))],
            ['Size:', SnpPhylogenyUtils.get_snp_matrix_size(snp_matrix)]
        ]
        section.add_table(table_data_snp_matrix, table_attributes=[('class', 'information')])

        # Add other output files table
        table_data = []
        for sample, output in output_files.items():
            row: list[Union[str, HtmlTableCell]] = [sample.name_full]
            for file_ in output:
                if file_ is None:
                    row.append('Not available')
                else:
                    relative_path = Path(sample.name_valid) / file_.name
                    row.append(HtmlTableCell('Download', link=str(section.add_file(file_, relative_path))))
            table_data.append(row)
        header = ['Sample'] + column_names
        section.add_table(table_data, header, [('class', 'data')])
        report.add_html_object(section)
        section.copy_files(report.output_dir)
        report.save()

    @staticmethod
    def add_metrics_section(report: HtmlReport, stats: dict[Sample, list], header: list[str]) -> None:
        """
        Adds a section with analysis metrics to the report.
        :param report: Report
        :param stats: Stats
        :param header: Table header
        :return: None
        """
        section = HtmlReportSection('Analysis metrics')

        # Add to output
        table_data = [[sample.name_full] + sts for sample, sts in stats.items()]
        section.add_table(table_data, header, [('class', 'data')])

        # Save as TSV
        data_metrics = pd.DataFrame([{k: v for k, v in zip(header, row)} for row in table_data])
        data_metrics.to_csv(Path(report.output_dir) / 'metrics.tsv', sep='\t', index=False)
        section.add_link_to_file('Download (TSV)', Path('metrics.tsv'))
        report.add_html_object(section)
        report.save()

    @staticmethod
    def get_snp_matrix_size(snp_matrix: Path) -> int:
        """
        Returns the length of the given SNP matrix.
        Note: this function only checks the length of the first entry.
        :param snp_matrix: SNP matrix
        :return: SNP matrix size
        """
        with snp_matrix.open() as handle:
            entries = list(SeqIO.parse(handle, 'fasta'))
        return len(entries[0].seq)

    @staticmethod
    def check_snp_matrix_size(snp_matrix: Path, size_min: int = 10, size_max: int = 25000) -> None:
        """
        Checks if the size of the SNP matrix is OK to perform model selection / tree building.
        :param snp_matrix: SNP matrix
        :param size_min: Minimum SNP matrix size
        :param size_max: Maximal SNP matrix size
        :return: None
        """
        snp_matrix_size = SnpPhylogenyUtils.get_snp_matrix_size(snp_matrix)
        if snp_matrix_size > size_max:
            raise ValueError(
                f'SNP matrix is too big ({snp_matrix_size}, max={size_max}) to perform model selection / tree building.')
        elif snp_matrix_size < size_min:
            raise ValueError(
                f'SNP matrix is too small ({snp_matrix_size}, min={size_min}) to perform model selection / tree building.')

    @staticmethod
    def add_model_selection_section(report: HtmlReport, model_selection: Optional[ModelSelection] = None,
                                    error_message: Optional[str] = None) -> None:
        """
        Adds the section with the model selection results.
        :param report: Report
        :param model_selection: Model selection tool
        :param error_message: Error message
        :return: None
        """
        section = HtmlReportSection('Model selection')
        if error_message is not None:
            section.add_error_message(error_message)
        else:
            relative_path = Path('model_selection_out.csv')
            section.add_file(model_selection.tool_outputs['CSV'][0].path, relative_path)
            table_data = [
                ['Selected model:', model_selection.informs['model_full']],
                ['Rates among sites:', model_selection.informs['rates_among_sites_full']],
                ['Overview:', HtmlTableCell('Download (CSV)', link=str(relative_path))]
            ]
            section.add_table(table_data, table_attributes=[('class', 'information')])
        report.add_html_object(section)
        section.copy_files(report.output_dir)
        report.save()

    @staticmethod
    def run_model_selection(snp_matrix: ToolIOFile, args: argparse.Namespace) -> ModelSelection:
        """
        Runs the MEGA model selection step.
        :param snp_matrix: SNP matrix
        :param args: Command line arguments
        :return: Model selection
        """
        model_selection = ModelSelection()
        model_selection.add_input_files({'FASTA': [snp_matrix]})
        MEGAUtils.update_model_selection_parameters(
            model_selection, args.missing_data, args.branch_swap, args.site_cov_cutoff, args.threads)
        working_dir = Path(args.working_dir) / 'model_selection'
        if not working_dir.exists():
            working_dir.mkdir(parents=True)
        model_selection.run(working_dir)
        return model_selection

    @staticmethod
    def add_tree_building_section(report: HtmlReport, newick_path: Optional[Path] = None,
                                  error_message: Optional[str] = None) -> None:
        """
        Adds the tree building section.
        :param report: Report
        :param newick_path: Path to the newick tree
        :param error_message: Error message
        :return: None
        """
        section = HtmlReportSection('Tree construction')
        if error_message is not None:
            section.add_error_message(error_message)
        else:
            # Add download link
            relative_path = Path('tree.nwk')
            section.add_file(newick_path, relative_path)
            table_data = [['Tree:', HtmlTableCell('Download (Newick)', link=str(relative_path))]]
            section.add_table(table_data, table_attributes=[('class', 'information')])

            # Render tree
            output_path = Path(report.output_dir) / 'tree.png'
            NewickUtils.create_image_figtree(
                newick_path, Path(str(files('camel').joinpath('resources/figtree/template_snp_tree.txt'))),
                output_path, 480, NewickUtils.calculate_tree_image_height(256, NewickUtils.count_leaves(newick_path))
            )
            section.add_html_object(HtmlElement('img', attributes=[('src', 'tree.png'), ('border', '1')]))

        report.add_html_object(section)
        section.copy_files(report.output_dir)
        report.save()

    @staticmethod
    def run_tree_building(snp_matrix: ToolIOFile, model: str, rates: str, args: argparse.Namespace) -> \
            MLTreeConstruction:
        """
        Builds a phylogenetic tree based on the given SNP matrix.
        :param snp_matrix: SNP matrix
        :param model: Selected model (e.g. 'T92')
        :param rates: Rates among sites (e.g. 'U')
        :param args: Command line arguments
        :return: Tree building tool instance
        """
        tree_building = MLTreeConstruction()
        tree_building.add_input_files({'FASTA': [snp_matrix]})
        MEGAUtils.update_tree_building_parameters(
            tree_building, model, rates, args.bootstraps, args.missing_data, args.site_cov_cutoff,
            args.ml_method, args.branch_swap, args.threads)
        working_dir = Path(args.working_dir) / 'tree_building'
        if not working_dir.exists():
            working_dir.mkdir(parents=True)
        tree_building.run(working_dir)
        return tree_building
