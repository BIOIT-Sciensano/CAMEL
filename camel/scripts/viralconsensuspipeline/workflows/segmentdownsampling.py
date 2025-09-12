import dataclasses
import json
import logging
from importlib.resources import files
from pathlib import Path
from typing import Union, Optional

from camel.app.camel import Camel
from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils


@dataclasses.dataclass
class SegmentDownsamplingOutput:
    """
    Holder for the output of the variant filtering output.
    """
    fq_out: FastqInput
    informs: list[dict]


class SegmentDownsamplingWorkflow:
    """
    Down samples the input dataset to a maximum coverage for each input segment.
    """

    def __init__(self, dir_: Path) -> None:
        """
        Initializes this workflow.
        :param dir_: Working directory
        :return: None
        """
        self._dir = dir_
        if not self._dir.exists():
            logging.info(f'Creating working directory: {self._dir}')
            self._dir.mkdir(parents=True)

    @staticmethod
    def calculate_ratio(depth: int, max_depth: int) -> Union[float, None]:
        """
        Calculates the ratio for downsampling.
        :param depth: Sample depth
        :param max_depth: Maximum depth
        :return: Ratio (None if no downsampling is needed)
        """
        if depth < max_depth:
            return None
        return max_depth / depth

    def run(self, bam_in: Path, input_type: str, json_mapping: Path, max_depth: int, bed_primers: Optional[Path] = None,
            threads: int = 8) -> SegmentDownsamplingOutput:
        """
        Runs the read mapping workflow.
        :param bam_in: Input BAM file
        :param input_type: Read type
        :param json_mapping: Mapping stats in JSON format
        :param max_depth: Maximum depth
        :param bed_primers: BED file with primer locations
        :param threads: Number of threads
        :return: Output holder
        """
        # Parse the mapping statistics
        with json_mapping.open() as handle:
            data_mapping = json.load(handle)

        # Create configfile
        path_config = SnakePipelineUtils.generate_config_file({
            'input': {'bam': str(bam_in), 'input_type': input_type},
            'downsampling': {
                seq_id: SegmentDownsamplingWorkflow.calculate_ratio(d['depth_median'], max_depth)
                for seq_id, d in data_mapping['by_chr'].items()},
            'bed_primers': str(bed_primers) if bed_primers is not None else None
        }, self._dir)
        path_snakefile = str(files('camel').joinpath(
            'scripts/viralconsensuspipeline/workflows/segmentdownsampling.smk'))

        # Run snakemake
        targets = {'fq': Path(f'merged/{input_type}/fq_dict.io'), 'informs': 'informs_all.json'}
        SnakePipelineUtils.run_snakemake(
            path_snakefile, str(path_config), list(targets.values()), working_dir=self._dir, threads=threads)

        # Collect output
        fq_dict_out = FastqInput.from_fq_dict(self._dir / targets['fq'], input_type)
        with (self._dir / targets['informs']).open() as handle:
            informs = json.load(handle)
        return SegmentDownsamplingOutput(fq_out=fq_dict_out, informs=informs)


if __name__ == '__main__':
    Camel.get_instance()
