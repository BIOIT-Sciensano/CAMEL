from pathlib import Path
from typing import Union

import vcf
# noinspection PyProtectedMember
from vcf.model import _Record as Record

from camel.app.loggers import logger
from camel.app.tools.variantfiltering.basefilter import BaseFilter


class DistanceFilter(BaseFilter):
    """
    Filters variants based on distance.

    Note:
        This code does not check for duplicate SNP positions. If this is the case the tool will only keep the last one
        in the VCF file.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Variant Filter: Distance', '0.1')

    @property
    def full_name(self) -> str:
        """
        Returns the full name for this filter.
        :return: Full name
        """
        return 'Distance'

    @property
    def description(self) -> str:
        """
        Returns the description for this filter.
        :return: Description
        """
        if 'keep_best' in self._parameters:
            extra = 'Variant with highest SNP quality is kept'
        else:
            extra = 'Both variants are removed'
        return 'Distance ≥<b>{}</b> between variants ({})'.format(
            self._parameters['min_distance'].value, extra)

    def _apply_filter(self) -> None:
        """
        Applies the filtering on the variants.
        :return: None
        """
        with open(self._tool_inputs['VCF_GZ'][0].path, 'rb') as handle:
            vcf_reader = vcf.Reader(handle)
            all_variants = list(vcf_reader)
        logger.info(f'Parsed {len(all_variants):,} variants')
        positions_to_filter = self.__get_positions_to_filter(vcf_reader, all_variants)
        logger.info(f'Filtered positions: {len(positions_to_filter):,}')
        bed_file = self.folder / 'positions_to_filter.bed'
        self.__create_bed_file(positions_to_filter, bed_file)
        logger.info('Created bed file of variants to filter.')
        self.__build_command(bed_file if len(positions_to_filter) > 0 else None)
        self._execute_command()

    def __build_command(self, bed_file: Union[Path, None]) -> None:
        """
        Builds the command for this tool.
        :param bed_file: Path to the BED file of the positions to soft-filter or filter.
        :return: None
        """
        file_parameter = '--targets-file ^' if 'soft_filter' not in self._parameters else '--mask-file '
        parts = [self._tool_command,
            str(self._tool_inputs['VCF_GZ'][0].path),
            '--output-type z',
            f'--output {self.output_path}'
        ]
        if bed_file is not None:
            parts.append(f"{file_parameter}{str(bed_file)}")
        self._command.command = ' '.join(parts + self._get_soft_filter_options())

    @staticmethod
    def __create_bed_file(positions_to_filter: list[tuple[str, int]], bed_file: Path) -> None:
        """
        Creates a BED file of variants to filter to pass to the bcftools filter command.
        :param positions_to_filter: The positions that need to be (soft-)filtered out.
        :param bed_file: Path to the BED file
        :return: None
        """
        # Sort positions_to_filter by chromosome and then by position
        positions_to_filter.sort()

        with bed_file.open('w') as handle:
            for chrom, pos in positions_to_filter:
                # Convert 1-based position to 0-based
                start_pos = pos - 1
                end_pos = pos  # BED format is 0-based, half-open interval
                handle.write(f"{chrom}\t{start_pos}\t{end_pos}\t.\t0\t+\n")

    def __get_positions_to_filter(self, vcf_reader: vcf.Reader, variants: list[Record]) -> list[tuple[str, int]]:
        """
        Returns a list with positions that should be filtered.
        :param vcf_reader: VCF reader
        :param variants: List of variants
        :return: List of positions that are removed by the filter
        """
        removed_positions = []
        interval_size = int(self._parameters['min_distance'].value)
        for contig_name, contig in vcf_reader.contigs.items():
            variants_by_pos = {v.POS: v for v in variants if v.CHROM == contig_name and (
                    v.FILTER is None) or len(v.FILTER) == 0}

            for i in range(0, contig.length - interval_size):
                interval = range(i, i + interval_size)
                variants_in_interval = [variants_by_pos[pos] for pos in interval if pos in variants_by_pos]

                # Do nothing if there is no SNP or only one
                if 0 <= len(variants_in_interval) < 2:
                    continue

                # Multiple SNPs - check which ones should be removed
                if 'keep_best' in self._parameters:
                    best_variant = max(variants_in_interval, key=lambda x: x.QUAL)
                    removed_variants = [v for v in variants_in_interval if v is not best_variant]
                else:
                    removed_variants = variants_in_interval
                removed_positions.extend([(v.CHROM, v.POS) for v in removed_variants])
        return list(set(removed_positions))
