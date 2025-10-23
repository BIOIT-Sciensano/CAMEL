import datetime
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Union, Any

import pandas as pd
import vcf
# noinspection PyProtectedMember
from vcf.model import _Record as VcfRecord

from camel.app.core.errors import InvalidToolInputError
from camel.app.loggers import logger
from camel.app.core.tool import Tool


@dataclass(frozen=True, unsafe_hash=True)
class Lineage:
    id_: str
    name: str
    main_spoligo: str
    rd_type: str

    @property
    def level(self) -> int:
        """
        Returns the level of the lineage, determined by the number of dots in the name.
        :return: Level
        """
        return self.id_.count('.')


@dataclass
class LineageSNP:
    chrom: str
    start: int
    end: int
    lineage: Lineage
    ref: str
    alt: str
    passes_filtering: bool = False


class SNPLineageDetector(Tool):
    """
    This tool is used to assign a SNP lineage to Mycobacterium tuberculosis complex strains.
    """

    def __init__(self):
        """
        Initializes this tool.
        """
        super().__init__('Mycobacterium: SNP lineage detector', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'VCF' not in self._tool_inputs:
            raise InvalidToolInputError('VCF input is required')
        if 'VCF_filt' not in self._tool_inputs:
            raise InvalidToolInputError("Filtered VCF input is required ('VCF_filt')")
        if 'BED' not in self._tool_inputs:
            raise InvalidToolInputError('BED input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        lineage_snps = self.__parse_lineage_snps()
        variants = SNPLineageDetector.__parse_vcf(self._tool_inputs['VCF'][0].path)
        variants_filtered = SNPLineageDetector.__parse_vcf(self._tool_inputs['VCF_filt'][0].path)
        logger.info(f"{len(variants_filtered)}/{len(variants)} (filtered) variants parsed.")

        # Check if variants are present in the database
        self._informs['detected_snps'] = self.__get_detected_snps(lineage_snps, variants, variants_filtered)
        counts_by_lineage = self.__count_lineage_snps()
        self._informs['counts_by_lineage'] = counts_by_lineage
        self._informs['detected_lineage_by_level'] = {i: self.__select_best_lineage_for_level(
            counts_by_lineage, i) for i in range(0, 5)}
        try:
            top_defined_lvl = max([k for k, l in self._informs['detected_lineage_by_level'].items() if l is not None])
            self._informs['detected_lineage'] = self._informs['detected_lineage_by_level'][top_defined_lvl]['lineage']
        except ValueError:
            self._informs['detected_lineage'] = Lineage('na', 'n/a', 'n/a', 'n/a')

        # Get the database version
        self._informs['db_version'] = self.__get_db_version()

    def __parse_lineage_snps(self) -> list[LineageSNP]:
        """
        Parses the file with the lineages.
        :return: Map of the positions to the corresponding lineage SNPs, Map of the lineages by name
        """
        data_lineage_snps = pd.read_table(self._tool_inputs['BED'][0].path, names=[
            'chr', 'start', 'end', 'lineage', 'ref', 'alt', 'lineage_name', 'spoligo', 'RD'])
        lineage_by_id = {}
        for id_, data in data_lineage_snps.fillna('n/a').groupby('lineage'):
            lineage_by_id[id_] = Lineage(
                id_=str(id_),
                name=data.iloc[0]['lineage_name'],
                main_spoligo=data.iloc[0]['spoligo'],
                rd_type=data.iloc[0]['RD'])
        logger.info(f'{len(lineage_by_id)} lineages parsed')

        lineage_snps = [LineageSNP(
            r['chr'], r['start'], r['end'], lineage_by_id[r['lineage']], r['ref'], r['alt'])
            for r in data_lineage_snps.to_dict('records')]
        logger.info(f"{len(lineage_snps)} lineage SNPs parsed")
        return lineage_snps

    @staticmethod
    def __parse_vcf(file_path: Path) -> list[VcfRecord]:
        """
        Parses the keys of the variants.
        :param file_path: VCF file path
        :return: Mapping of variants per position
        """
        with open(file_path) as handle:
            vcf_reader = vcf.Reader(handle)
            return list(vcf_reader)

    def __get_detected_snps(self, snps: list[LineageSNP], variants: list[VcfRecord], variants_filt: list[VcfRecord]) \
            -> list[LineageSNP]:
        """
        Returns the SNP positions detected in the sample.
        :param snps: Lineage SNP positions
        :param variants: Variants
        :param variants_filt: Filtered variants
        :return: Lineage SNPs present in the sample
        """
        variants_by_pos = {v.POS: v for v in variants}
        variants_filt_by_pos = {v.POS: v for v in variants_filt}
        detected_snps = []
        for snp in snps:
            if snp.lineage.id_ in ('lineage4', 'lineage4.9'):
                # Position needs to match reference
                if snp.start not in variants_by_pos:
                    detected_snps.append(snp)
            else:
                # Position needs to match DB snp
                if (snp.start in variants_by_pos) and (str(variants_by_pos[snp.start].ALT[0]) == snp.alt):
                    detected_snps.append(snp)
                    snp.passes_filtering = snp.start in variants_filt_by_pos
        return detected_snps

    def __count_lineage_snps(self) -> dict[Lineage, int]:
        """
        Counts the number of SNPs for each lineage.
        :return: Dictionary of number of SNPs for each lineage
        """
        counts_by_lineage = {}
        for snp in self._informs['detected_snps']:
            if snp.lineage not in counts_by_lineage:
                counts_by_lineage[snp.lineage] = 0
            counts_by_lineage[snp.lineage] += 1
        return counts_by_lineage

    def __select_best_lineage_for_level(self, count_by_lineages: dict[Lineage, int], level: int) \
            -> Union[dict[Any, int], None]:
        """
        Selects the best lineage for the given level.
        :param count_by_lineages: Count by lineage
        :param level: Lineage level
        :return: Lineage (if there is one at the given level)
        """
        filtered_counts = [{'lineage': l, 'count': c} for l, c in count_by_lineages.items() if l.level == level]
        if len(filtered_counts) == 0:
            return None
        return sorted(filtered_counts, key=lambda x: x['count'], reverse=True)[0]

    def __get_db_version(self) -> str:
        """
        Extracts the database version.
        :return: Database version
        """
        m = re.match(r'.*-(\d{4})_(\d{2})_(\d{2})\.bed', self._tool_inputs['BED'][0].path.name)
        if not m:
            logger.warning(f'Cannot determine database version')
            return 'n/a'
        date = datetime.date(year=int(m.group(1)), month=int(m.group(2)), day=int(m.group(3)))
        return date.strftime('%d/%m/%Y')
