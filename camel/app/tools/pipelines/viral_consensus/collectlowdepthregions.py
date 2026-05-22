from pathlib import Path

import pandas as pd
from Bio import SeqIO
from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.tool import Tool
from camel.app.loggers import logger


class CollectLowDepthRegions(Tool):
    """
    Collects low-depth regions from an input BAM file.
    It uses bedtools to determine the coverage based on the input BAM file.
    Then it uses bedtools complement to determine low-depth regions.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Collect low depth regions', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA input is required')
        if 'BAM' not in self._tool_inputs:
            raise InvalidToolInputError('BAM input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        """
        path_bed_gcov = self._run_bedtools_genome_cov(self._tool_inputs['BAM'][0].path)
        path_bed_low_cov = self._extract_low_coverage_regions(path_bed_gcov)
        path_bed_merged = self._merge_intervals(path_bed_low_cov)
        if int(self._parameters['gap_len_cutoff'].value) > 0:
            path_bed_merged = self._remove_small_intervals(
                path_bed_merged, int(self._parameters['gap_len_cutoff'].value))
        self._extract_gap_stats(path_bed_merged)
        self._tool_outputs['BED'] = [ToolIOFile(path_bed_merged)]

    def _run_bedtools_genome_cov(self, path_bam: Path) -> Path:
        """
        Runs bedtools genomecov to determine per position depth.
        :param path_bam: Input BAM file
        :return: Path to output file
        """
        path_bed_gcov = self.folder / 'gcov.bed'
        self._command.command = ' '.join([
            self._build_dependencies(), 'bedtools genomecov', f"-ibam {path_bam}", '-bga', f'> {path_bed_gcov}'])
        self._execute_command()
        if not self._command.exit_code == 0:
            raise ToolExecutionError(f'Error running bedtools genome cov: {self._command.stderr}')
        logger.info(f'bedtools genome cov output file generated: {path_bed_gcov}')
        return path_bed_gcov

    def _extract_low_coverage_regions(self, path_bed: Path) -> Path:
        """
        Extracts the low coverage regions.
        :param path_bed: Input BED file with depth statistics
        :return: Path to BED file with low coverage regions
        """
        data_depth = pd.read_table(
            path_bed, names=['chr', 'start', 'stop', 'depth'], keep_default_na=False, na_values='-')
        logger.info(f'Depth data parsed ({len(data_depth):,} rows)')
        is_low_cov = data_depth['depth'] < int(self._parameters['gap_depth_cutoff'].value)
        data_depth_low_cov = data_depth[is_low_cov].copy()
        data_depth_low_cov.sort_values(by=['chr', 'start'], inplace=True)
        path_bed_low_cov = self.folder / 'low_cov.bed'
        data_depth_low_cov.to_csv(path_bed_low_cov, sep='\t', index=False, header=False)
        logger.info(f'{len(data_depth_low_cov):,}/{len(data_depth):,} rows below depth threshold')
        logger.info(f'BED file with low depth regions created: {path_bed_low_cov}')
        return path_bed_low_cov

    def _merge_intervals(self, path_bed: Path) -> Path:
        """
        Merges the intervals in the input BED file.
        :param path_bed: Input BED file
        :return: Path to BED file with merged intervals
        """
        path_bed_merged = self.folder / 'low_cov_merged.bed'
        command = Command(' '.join([
            self._build_dependencies(), f'bedtools merge -i {path_bed} > {path_bed_merged}']))
        command.run(self.folder)
        if not command.exit_code == 0:
            raise ToolExecutionError(f'Error executing bedtools: {command.stderr}')
        logger.info(f'BED file with merged regions created: {path_bed_merged}')
        return path_bed_merged

    def _remove_small_intervals(self, path_bed: Path, min_size: int) -> Path:
        """
        Removes the small intervals (i.e., small deletions) from the input BED file.
        :param path_bed: Input BED file
        :param min_size: Minimum interval size
        :return: Updated BED file
        """
        data_intervals = pd.read_table(path_bed, names=['chr', 'start', 'end'], keep_default_na=False, na_values='-')
        interval_size = data_intervals['end'] - data_intervals['start']
        data_intervals_filt = data_intervals[interval_size >= min_size]
        logger.info(f'{len(data_intervals_filt)}/{len(data_intervals)} intervals above size limit (>={min_size})')
        path_out = self.folder / 'low_cov_merged_filt.bed'
        data_intervals_filt.to_csv(path_out, sep='\t', index=False, header=False)
        return path_out

    def _extract_gap_stats(self, path_bed_merged: Path) -> None:
        """
        Extracts the statistics for the detected gaps.
        :param path_bed_merged: Merged BED file
        """
        # Get FASTA lengths
        with open(self._tool_inputs['FASTA'][0].path) as handle:
            len_by_seq_id = {s.id: len(s) for s in SeqIO.parse(handle, 'fasta')}

        # Parse BED file
        data_depth_merged = pd.read_table(
            path_bed_merged, names=['chr', 'start', 'end'], keep_default_na=False, na_values='-')
        data_depth_merged['size'] = data_depth_merged['end'] - data_depth_merged['start']

        # Collect statistics (Global)
        self._informs['total'] = {
            'nb_gaps': len(data_depth_merged),
            'total_gap_size': sum(data_depth_merged['size']),
            'length': sum(len_by_seq_id.values()),
            'perc_gaps': 100 * sum(data_depth_merged['size']) / sum(len_by_seq_id.values())
        }

        # Collect statistics (by sequence id)
        self._informs['by_seq_id'] = {}
        data_depth_by_seq_id = {seq_id: data for seq_id, data in data_depth_merged.groupby('chr')}
        for seq_id, length in len_by_seq_id.items():
            self._informs['by_seq_id'][seq_id] = {
                'nb_gaps': len(data_depth_by_seq_id.get(seq_id, [])),
                'total_gap_size': sum(data_depth_by_seq_id[seq_id]['size']) if seq_id in data_depth_by_seq_id else 0,
                'perc_gaps': 100 * (sum(
                    data_depth_by_seq_id[seq_id]['size']) if seq_id in data_depth_by_seq_id else 0) / length,
                'length': length
            }
