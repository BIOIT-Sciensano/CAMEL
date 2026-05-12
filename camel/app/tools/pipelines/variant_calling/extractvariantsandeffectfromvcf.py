from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from vcf.model import _Record as VcfRecord

from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool
from camel.app.core.utils import toolutils
from camel.app.core.utils.vcfutils import retrieve_variants
from camel.app.loggers import logger
from camel.app.toolkits.export.tsvexporter import TsvExporter


@dataclass
class Variant:
    """
    This class represents a variant from a VCF file.
    """
    position: int
    type_of_variant: str
    variant: str
    effect: str
    gene: str
    allele_freq: float
    quality: float

    def to_list(self) -> list:
        """
        Returns the variant information as a list.
        :return: list
        """
        return [self.position, self.type_of_variant, self.variant, self.effect, self.gene, self.allele_freq,
                self.quality]


class ExtractVariantsAndEffectFromVCF(Tool):
    """
    This tool is used to extract variants from a VCF file and format them in a pandas table.
    The output of the table consists of the following columns:
    Position, Type of variant, the genomic change of the variant (e.g., which mutation was identified),
    Effect (if BCSQ field is present), Gene (if BCSQ field is present), Allele Frequency, QUAL score
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__("Extract Variants From VCF", "0.1")

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['VCF'])

    def _execute_tool(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        variants_dataframe = self.__parse_variants_for_output_table(self._tool_inputs['VCF'][0].path)
        table_path = self._folder / 'all_variants.tsv'
        TsvExporter.export(variants_dataframe.values.tolist(), list(variants_dataframe.columns), table_path)
        self._tool_outputs['TSV'] = [ToolIOFile(table_path)]

    @staticmethod
    def __parse_effect(vcf_record: VcfRecord) -> tuple[str | None, str | None]:
        """
        Parses the mutation effect from the CSQ annotation.
        Note: only extracts it for protein coding regions
        :param vcf_record: Input record
        :return: Mutation effect and Gene
        """
        # Check if BCSQ annotation is present
        if 'BCSQ' not in vcf_record.INFO:
            logger.warning(f'BCSQ info missing for: {vcf_record.CHROM}:{vcf_record.POS}')
            return None, None

        # Parse annotation
        parts = vcf_record.INFO['BCSQ'][0].split('|')
        if parts[0].startswith('&'):
            return None, None
        if parts[0].startswith('@'):
            return parts[0], '-'
        return parts[0], parts[1]

    def __parse_variants_for_output_table(self, input_vcf: Path) -> pd.DataFrame:
        """
        Parses the variants list for the summary variant table in the report.
        :param input_vcf: Input VCF file
        :return: Pandas DF of all variants with their associated effect
        """
        output_dictionary = {}
        positions_to_check_at_the_end = {}
        var_list = retrieve_variants(input_vcf)
        for var in var_list:
            effect, gene = self.__parse_effect(var)
            variant = f'{var.REF}->{var.ALT[0]}'
            type_of_var = 'Indel' if var.INFO.get('INDEL', False) is True else var.var_type
            if effect is None:
                effect, gene = 'Unknown', 'Unknown'
            output_dictionary[var.POS] = Variant(position=var.POS, type_of_variant=type_of_var,
                                                 variant=variant, effect=effect, gene=gene,
                                                 allele_freq=var.INFO.get('AF', 0), quality=var.QUAL)

            # This condition checks whether the variant effect is a copy of another (symbolized by the '@')
            # I'm therefore storing the position the VCF entry refers to
            if effect.startswith('@'):
                positions_to_check_at_the_end[var.POS] = int(effect[1:])

        # For each position with an @, I replace the empty effect and gene with the referenced ones
        for unknown_pos, ref_pos in positions_to_check_at_the_end.items():
            output_dictionary[unknown_pos].effect = output_dictionary[ref_pos].effect
            output_dictionary[unknown_pos].gene = output_dictionary[ref_pos].gene
        output_dataframe = pd.DataFrame(data=[var.to_list() for key, var in output_dictionary.items()],
                                        columns=['Position', 'Type', 'Variant', 'Effect', 'Gene', 'AF', 'Quality'])
        return output_dataframe
