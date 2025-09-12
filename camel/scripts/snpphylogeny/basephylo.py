import abc
import argparse
from collections.abc import Sequence
from typing import Optional

from camel.app.components.phylogeny.snpphylogenyutils import (
    MappingInput,
    Sample,
    SnpPhylogenyUtils,
)
from camel.app.error import InvalidToolInputError
from camel.app.error import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.mega.modelselection import ModelSelection
import sys


class BasePhylo(metaclass=abc.ABCMeta):
    """
    Base class for the SNP phylogeny pipelines.
    """

    def __init__(self, pipeline_name: str, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the base class.
        :param pipeline_name: Pipeline name (e.g. 'Samtools')
        :param args: Command line arguments
        """
        self._pipeline_name = pipeline_name
        self._args = self._parse_arguments(args)
        self._report = SnpPhylogenyUtils.initialize_report(self._pipeline_name, self._args)
        self._samples = self.__extract_samples()
        self._informs = []

    @property
    def samples_by_name(self) -> dict[str, Sample]:
        """
        Returns the samples as a dictionary with the sample name as key.
        :return: Samples by name
        """
        return {s.name_valid: s for s in self._samples}

    @staticmethod
    @abc.abstractmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        pass

    def __extract_samples(self) -> list[Sample]:
        """
        Extracts sample objects from the provided input.
        :return: List of samples
        """
        try:
            return SnpPhylogenyUtils.extract_samples(self._args)
        except InvalidToolInputError as err:
            logger.error(f"Invalid input: {err}")
            self._report.add_error_message(str(err))
            self._report.save()
            exit(0)

    def _get_mapping_input(self) -> dict[Sample, MappingInput]:
        """
        Returns the input for the read mapping.
        :return: Mapping input per sample
        """
        fq_by_sample = SnpPhylogenyUtils.symlink_input_files(self._samples, self._args.working_dir)
        if self._args.trim_reads:
            trimming_output_by_sample = SnpPhylogenyUtils.trim_all_reads(
                fq_by_sample, self._args.working_dir / 'trimming', self._args.adapter, self._args.threads)
            SnpPhylogenyUtils.add_trimming_section(self._report, trimming_output_by_sample)
            mapping_input_by_sample = {}
            for sample, output in trimming_output_by_sample.items():
                mapping_input_by_sample[sample] = MappingInput(
                    pe=output.trimmed_reads_pe,
                    se_fwd=output.trimmed_reads_se_fwd[0] if len(output.trimmed_reads_se_fwd) > 0 else None,
                    se_rev=output.trimmed_reads_se_rev[0] if len(output.trimmed_reads_se_rev) > 0 else None
                )
            self._informs.append(trimming_output_by_sample[self._samples[0]].informs_trimming)
            return mapping_input_by_sample
        else:
            SnpPhylogenyUtils.add_trimming_section_empty(self._report)
            return {s: MappingInput(pe=fq) for s, fq in fq_by_sample.items()}

    def _run_model_selection(self, snp_matrix: ToolIOFile) -> ModelSelection:
        """
        Runs the model selection.
        Quits the programs if the SNP matrix is too small.
        :param snp_matrix: SNP matrix
        :return: Model selection tool instance
        """
        try:
            SnpPhylogenyUtils.check_snp_matrix_size(snp_matrix.path)
        except ValueError as err:
            SnpPhylogenyUtils.add_model_selection_section(self._report, error_message=str(err))
            sys.exit(0)
        else:
            model_selection = SnpPhylogenyUtils.run_model_selection(snp_matrix, self._args)
            SnpPhylogenyUtils.add_model_selection_section(self._report, model_selection=model_selection)
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
            tree_building = SnpPhylogenyUtils.run_tree_building(
                snp_matrix, model_selection.informs['model'], model_selection.informs['rates_among_sites'], self._args)
            SnpPhylogenyUtils.add_tree_building_section(self._report, tree_building.tool_outputs['NWK'][0].path)
            self._informs.append(tree_building.informs)
        except ToolExecutionError:
            SnpPhylogenyUtils.add_tree_building_section(
                self._report, error_message='Error constructing tree, SNP matrix might be too small')
