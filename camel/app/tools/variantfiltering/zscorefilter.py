import math
from pathlib import Path

from camelcore.app.command import Command
from cyvcf2 import VCF, Variant, Writer

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
        y_mult_val = self.get_param_value('y_multiplier')
        y_mult = y_mult_val if y_mult_val is not None else 'n/a'
        return f'Z-score ≥<b>{self._parameters["min_zscore"].value}</b> and Y-multiplier ≥<b>{y_mult}</b>...'

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
        denominator = math.sqrt(x + y) if (x + y) > 0 else 1
        return float(x - y) / denominator

    @staticmethod
    def get_actg_counts(pileup_line: str) -> list[int]:
        """
        Get the count for each base
        :param pileup_line: Line of the pileup output.
        :return: A, C, T, G counts
        """
        return [pileup_line.upper().count(base) for base in ('A', 'C', 'T', 'G')]

    def __get_actg_counts_by_position(
        self, variants: list[Variant]
    ) -> dict[tuple[str, int], list[int]]:
        """
        Creates a pileup file with every position that covers a variant.
        :param variants: List of variants
        :return: Path to pileup file
        """
        positions_file = Path(self._folder) / 'positions.txt'
        with open(positions_file, 'w') as handle_out:
            for variant in variants:
                handle_out.write(f'{variant.CHROM}\t{variant.POS}\n')

        path_pileup = Path(self._folder) / 'positions.pileup'
        pileup_command = Command(' '.join([
            self._build_dependencies(),
            'samtools mpileup',
            '--count-orphans',
            '--positions', str(positions_file),
            str(self._tool_inputs['BAM'][0].path),
            f'> {path_pileup}'
        ]))
        pileup_command.run(self._folder)

        # Parse the output
        actg_by_pos = {}
        with path_pileup.open() as handle:
            for line in handle.readlines():
                parts = line.strip().split('\t')
                actg_by_pos[(parts[0], int(parts[1]))] = ZScoreFilter.get_actg_counts(
                    parts[4]
                )
        return actg_by_pos

    def __get_pos_to_filter(
        self, actg_by_pos: dict[tuple[str, int], list[int]], all_variants: list[Variant]
    ) -> set[tuple[str, int]]:
        """
        Filters all variant positions based on the Z-score.
        :param actg_by_pos: ACTG counts by position
        :param all_variants: List of all variants
        :return: List with the positions that should be removed.
        """
        filtered_positions = []
        for variant in all_variants:
            if variant.is_indel or (
                variant.FILTER is not None and len(variant.FILTER) > 0
            ):
                continue
            actg_counts = actg_by_pos[(variant.CHROM, variant.POS)]
            x = max(actg_counts)
            y = sum(actg_counts) - x
            z_score = ZScoreFilter.calculate_zscore(x, y)

            if z_score < float(self._parameters['min_zscore'].value):
                filtered_positions.append((variant.CHROM, variant.POS))
                continue

            if (
                'y_multiplier' in self._parameters
                and x < float(self.get_param_value('y_multiplier')) * y
            ):
                filtered_positions.append((variant.CHROM, variant.POS))
                continue

        return set(filtered_positions)

    def _apply_filter(self) -> None:
        """
        Applies the filtering on the variants.
        :return: None
        """
        with VCF(str(self._tool_inputs['VCF_GZ'][0].path)) as vcf_reader:
            all_variants: list[Variant] = list(vcf_reader)
        actg_counts_by_position = self.__get_actg_counts_by_position(all_variants)
        filtered_positions = self.__get_pos_to_filter(
            actg_counts_by_position, all_variants
        )
        output_uncompressed = self.__create_output_file(
            self._tool_inputs['VCF_GZ'][0].path, filtered_positions
        )
        command_bgzip = Command(f'bgzip {output_uncompressed}')
        self._execute_command(command_bgzip)

    def __create_output_file(
        self,
        path_in: Path,
        positions_failing_filt: set[tuple[str, int]],
    ) -> Path:
        """
        Creates the output file.
        :return: Output file (uncompressed VCF)
        """
        vcf_reader = VCF(str(path_in))

        # Add custom filter
        if self.get_param_value('soft_filter') is True:
            vcf_reader.add_filter_to_header({
                'ID': 'z_score',
                'Description': 'Z-score as described in CSI phylogeny (custom)'
            })

        # Write updated VCF file
        output_uncompressed = self.output_path.parent / self.output_path.name.replace('.gz', '')
        writer = Writer(str(output_uncompressed), vcf_reader)
        for variant in vcf_reader:
            if variant.is_snp and ((variant.CHROM, variant.POS) in positions_failing_filt):
                if self.get_param_value('soft_filter') is True:
                    if variant.FILTER in (None, '.', 'PASS'):
                        variant.FILTER = 'z_score'
                    else:
                        variant.FILTER = f"{variant.FILTER};z_score"
                else:
                    continue
            writer.write_record(variant)

        writer.close()
        vcf_reader.close()
        return output_uncompressed
