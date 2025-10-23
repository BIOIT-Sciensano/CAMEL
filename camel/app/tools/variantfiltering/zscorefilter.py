import math
from pathlib import Path

import vcf
# noinspection PyProtectedMember
from vcf.model import _Record as Record
# noinspection PyProtectedMember
from vcf.parser import _Filter as Filter

from camel.app.core.command import Command
from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.variantfiltering.basefilter import BaseFilter


class ZScoreFilter(BaseFilter):
    """
    Filters variants based on Z-score.
    The Z-score is calculated as:
    Z = (X-Y) / sqrt(X+Y)
    Where X is the number of reads having the most common nucleotide at that position and Y the number of reads
    supporting other nucleotides.

    If the 'y_multiplier' parameter is set, positions for which the following condition does not hold are filtered out:
    X = Y * {y_multiplier}
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Variant Filter: Z-score', '0.1')

    @property
    def full_name(self) -> str:
        """
        Returns the full name for this filter.
        :return: Full name
        """
        return 'Z-score'

    @property
    def description(self) -> str:
        """
        Returns the description for this filter.
        :return: Description
        """
        return 'Z-score ≥<b>{}</b> and Y-multiplier ≥<b>{}</b> at variant position (see citation Kaas et al.)'.format(
            self._parameters['min_zscore'].value, self._parameters['y_multiplier'].value)

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if 'BAM' not in self._tool_inputs:
            raise InvalidToolInputError("No BAM input found")
        super()._check_input()

    @staticmethod
    def calculate_zscore(x: int, y: int) -> float:
        """
        Calculates the Z-score.
        :param x: Number of reads with the most common nucleotide.
        :param y: Number of other reads
        :return: Z-score
        """
        return float(x - y) / math.sqrt(x + y)

    @staticmethod
    def get_actg_counts(pileup_line: str) -> list[int]:
        """
        Get the count for each base
        :param pileup_line: Line of the pileup output.
        :return: A, C, T, G counts
        """
        return [pileup_line.upper().count(base) for base in ('A', 'C', 'T', 'G')]

    def __get_actg_counts_by_position(self) -> dict[tuple[str, int], list[int]]:
        """
        Creates a pileup file with every position that covers a variant.
        :return: Path to pileup file
        """
        positions_file = Path(self._folder) / 'positions.txt'
        with open(self._tool_inputs['VCF_GZ'][0].path, 'rb') as handle_in:
            with open(positions_file, 'w') as handle_out:
                for variant in vcf.Reader(handle_in):
                    handle_out.write(f'{variant.CHROM}\t{variant.POS}\n')

        path = Path(self._folder) / 'positions.pileup'
        pileup_command = Command('{}samtools mpileup --count-orphans --positions {} {} > {}'.format(
            self._build_dependencies(), positions_file, self._tool_inputs['BAM'][0].path, path))
        pileup_command.run(self._folder)

        actg_by_pos = {}
        with path.open() as handle:
            for line in handle.readlines():
                parts = line.strip().split('\t')
                actg_by_pos[(parts[0], int(parts[1]))] = ZScoreFilter.get_actg_counts(parts[4])
        return actg_by_pos

    def __get_filtered_positions(self, actg_by_pos: dict[tuple[str, int], list[int]],
                                 all_variants: list[Record]) -> list[tuple[str, int]]:
        """
        Filters all variant positions based on the Z-score.
        :param actg_by_pos: ACTG counts by position
        :param all_variants: List of all variants
        :return: List with the positions that are kept.
        """
        filtered_positions = []
        for variant in all_variants:
            if variant.is_indel or (variant.FILTER is not None and len(variant.FILTER) > 0):
                continue
            actg_counts = actg_by_pos[(variant.CHROM, variant.POS)]
            x = max(actg_counts)
            actg_counts.remove(x)
            y = sum(actg_counts)
            z_score = float(x - y) / max(math.sqrt(x + y), 1)

            if z_score < float(self._parameters['min_zscore'].value):
                filtered_positions.append((variant.CHROM, variant.POS))
                continue

            if 'y_multiplier' in self._parameters and x < float(self._parameters['y_multiplier'].value) * y:
                filtered_positions.append((variant.CHROM, variant.POS))
                continue

        return filtered_positions

    def _apply_filter(self) -> None:
        """
        Applies the filtering on the variants.
        :return: None
        """
        with open(self._tool_inputs['VCF_GZ'][0].path, 'rb') as handle:
            vcf_reader = vcf.Reader(handle)
            all_variants = list(vcf_reader)
        actg_counts_by_position = self.__get_actg_counts_by_position()
        filtered_positions = self.__get_filtered_positions(actg_counts_by_position, all_variants)
        output_uncompressed = self.__create_output_file(vcf_reader, all_variants, filtered_positions)
        self._command = Command(f'bgzip {output_uncompressed}')
        self._execute_command()

    def __create_output_file(self, vcf_reader: vcf.Reader, all_variants: list[Record],
                             filtered_positions: list[tuple[str, int]]) -> Path:
        """
        Creates the output file.
        :return: Output file (uncompressed VCF)
        """
        vcf_reader.filters['z_score'] = Filter('z_score', 'Z-score as described in CSI phylogeny (custom)')
        output_uncompressed = self.output_path.parent / self.output_path.name.replace('.gz', '')
        with output_uncompressed.open('w') as handle:
            writer = vcf.Writer(handle, vcf_reader)
            for variant in all_variants:
                if (variant.CHROM, variant.POS) in filtered_positions:
                    if 'soft_filter' in self._parameters:
                        variant.FILTER = 'z_score'
                    else:
                        continue
                writer.write_record(variant)
            writer.close()
        return output_uncompressed
