import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake import snakemakeutils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.bcftools.bcftoolsview import BcftoolsView
from camel.resources.snakefile import variant_filtering, variant_calling


@dataclass
class VariantFilteringOutput:
    """
    This dataclass holds the output of the variant filtering workflow.
    """
    vcf_filtered: ToolIOFile
    stats: dict
    nb_of_variants: int
    informs: list[dict[str, Any]]


class VariantFilteringWrapper:
    """
    This class is used as a wrapper class around the variant filtering Snakemake workflow.
    """

    def __init__(self, working_dir: Path) -> None:
        """
        Initializes the variant calling wrapper.
        :param working_dir: Working directory
        """
        self._working_dir = working_dir
        self._output = None

    @property
    def output(self) -> VariantFilteringOutput:
        """
        Returns the filtered VCF file
        :return: VCF file path
        """
        return self._output

    def __convert_to_vcf_gz(self, vcf_file: Path) -> Path:
        """
        Converts an input VCF file to an indexed VCF_GZ file.
        :param vcf_file: Input VCF file
        :return: Indexed VCF_GZ file
        """
        bcftools_view = BcftoolsView()
        bcftools_view.add_input_files({'VCF': [ToolIOFile(vcf_file)]})
        input_dir = Path(self._working_dir, 'input')
        if not input_dir.is_dir():
            input_dir.mkdir(parents=True)
        bcftools_view.update_parameters(output_type='z')
        bcftools_view.run(input_dir)
        return bcftools_view.tool_outputs['VCF_GZ'][0].path

    def __create_input(self, vcf_gz_file: Path, bam_file: Path, input_type: str) -> None:
        """
        Creates the input files for the workflow.
        :param vcf_gz_file: Input VCF GZ file
        :param bam_file: Input BAM file
        :param input_type: Input type
        :return: None
        """
        for path, destination in [
            (vcf_gz_file, variant_calling.OUTPUT_UNFILTERED_VCF_GZ),
            (bam_file, variant_calling.get_bam({'working_dir': self._working_dir, 'input_type': input_type}))]:
            target_file = Path(self._working_dir, destination)
            if not target_file.parent.exists():
                target_file.parent.mkdir(parents=True)
            snakemakeutils.dump_object([ToolIOFile(path)] if path is not None else [], target_file)

    def run_workflow(
            self, sample_name: str, vcf_file: Path, bam_file: Path, input_type: str = 'illumina',
            filtering_options: dict = None, cores: int = 8) -> None:
        """
        Runs the variant calling workflow.
        :param sample_name: Sample name
        :param vcf_file: Input VCF file
        :param bam_file: Input BAM file
        :param input_type: Input type ('illumina' or 'ont')
        :param cores: Number of cores
        :param filtering_options: Dict
        :return: None
        """
        if not self._working_dir.exists():
            self._working_dir.mkdir(parents=True)
        vcf_gz_file = self.__convert_to_vcf_gz(vcf_file)
        self.__create_input(vcf_gz_file, bam_file, input_type)

        # Create config
        config_path = SnakePipelineUtils.generate_config_file({
            'working_dir': str(self._working_dir),
            'variant_filtering': filtering_options,
            'input_type': input_type,
            'sample_name': sample_name
        }, self._working_dir)

        # Execute Snakemake
        output_files = {
            'VCF': variant_filtering.OUTPUT_VCF,
            'STATS': variant_filtering.OUTPUT_STATS,
            'INFORMS': variant_filtering.OUTPUT_INFORMS_ALL
        }
        SnakePipelineUtils.run_snakemake(
            variant_filtering.SNAKEFILE, config_path, [Path(x) for x in output_files.values()], self._working_dir, cores)
        self.__set_output(output_files)

    def __set_output(self, output_files: dict[str, Path]) -> None:
        """
        Collects the output of the workflow.
        :param output_files: Output files
        :return: None
        """
        json_file = snakemakeutils.load_object(self._working_dir / output_files['STATS'])[0].path
        with open(json_file) as handle:
            stats = json.load(handle)
        self._output = VariantFilteringOutput(
            vcf_filtered=snakemakeutils.load_object(self._working_dir / output_files['VCF'])[0],
            stats=stats,
            nb_of_variants=stats['zscore']['variants_out'],
            informs=snakemakeutils.load_object(self._working_dir / output_files['INFORMS'])
        )
