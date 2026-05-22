from pathlib import Path

from cyvcf2 import VCF, Variant

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
        return f'Distance ≥<b>{self.get_param_value("min_distance")}</b> between variants ({extra})'

    def _apply_filter(self) -> None:
        """
        Applies the filtering on the variants.
        :return: None
        """
        with VCF(str(self._tool_inputs['VCF_GZ'][0].path)) as vcf_reader:
            all_variants: list[Variant] = list(vcf_reader)
            length_by_contigs = {
                c: l for c, l in zip(vcf_reader.seqnames, vcf_reader.seqlens)
            }
        logger.info(f'Parsed {len(all_variants):,} variants')
        positions_to_filter = self.__get_positions_to_filter(
            length_by_contigs, all_variants
        )
        logger.info(f'Filtered positions: {len(positions_to_filter):,}')
        bed_file = self.folder / 'positions_to_filter.bed'
        self.__create_bed_file(positions_to_filter, bed_file)
        logger.info('Created bed file of variants to filter.')
        self.__build_command(bed_file if len(positions_to_filter) > 0 else None)
        self._execute_command()

    def __build_command(self, bed_file: Path | None) -> None:
        """
        Builds the command for this tool.
        :param bed_file: Path to the BED file of the positions to soft-filter or filter.
        :return: None
        """
        file_parameter = (
            '--targets-file ^'
            if 'soft_filter' not in self._parameters
            else '--mask-file '
        )
        parts = [
            self._tool_command,
            str(self._tool_inputs['VCF_GZ'][0].path),
            '--output-type z',
            f'--output {self.output_path}',
        ]
        if bed_file is not None:
            parts.append(f"{file_parameter}{str(bed_file)}")
        self._command.command = ' '.join(parts + self._get_soft_filter_options())

    @staticmethod
    def __create_bed_file(
        positions_to_filter: list[tuple[str, int]], bed_file: Path
    ) -> None:
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

    def __get_positions_to_filter(
            self, len_by_contig: dict[str, int], variants: list[Variant]
    ) -> list[tuple[str, int]]:
        """
        Returns a list with positions that should be filtered.
        :param len_by_contig: Length of each contig
        :param variants: List of variants
        :return: List of positions that are removed by the filter
        """
        removed_positions = []
        interval_size = int(self._parameters['min_distance'].value)

        for contig_name, contig_len in len_by_contig.items():
            # Filter variants for this contig
            contig_variants = [
                v for v in variants
                if v.CHROM == contig_name and (v.FILTER is None or len(v.FILTER) == 0)
            ]

            # Sort by position for efficiency
            contig_variants.sort(key=lambda x: x.POS)

            # Process only adjacent variants
            for j, variant in enumerate(contig_variants):
                for k in range(j + 1, len(contig_variants)):
                    next_variant = contig_variants[k]
                    distance = next_variant.POS - variant.POS

                    if distance >= interval_size:
                        break  # No more nearby variants

                    # Found close variants
                    if 'keep_best' in self._parameters:
                        if variant.QUAL < next_variant.QUAL:
                            removed_positions.append((variant.CHROM, variant.POS))
                        else:
                            removed_positions.append((next_variant.CHROM, next_variant.POS))
                    else:
                        removed_positions.extend([
                            (variant.CHROM, variant.POS),
                            (next_variant.CHROM, next_variant.POS)
                        ])

        return list(set(removed_positions))
