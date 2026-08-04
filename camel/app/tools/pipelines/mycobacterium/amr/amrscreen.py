import collections
import json
import re
from pathlib import Path
from typing import Union

import pandas as pd
from Bio.SeqUtils import seq3
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import vcfutils
from cyvcf2 import VCF, Variant

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool
from camel.app.loggers import logger


class AMRScreen(Tool):
    """
    Screens the mutations from a VCF file against a database with mutations.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('AMR mutation screen', '0.1')
        self._variant_by_key: dict[tuple[int, str, str], str] | None = None
        self._data_regions: pd.DataFrame | None = None
        self._ab_short_by_name: dict[str, str] | None = None
        self._association_by_variant = None
        self._data_antibiotics: pd.DataFrame | None = None

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'VCF' not in self._tool_inputs:
            raise InvalidToolInputError("Mutation input is required (VCF)")
        if 'VCF_filt' not in self._tool_inputs:
            raise InvalidToolInputError("Filtered mutation input is required (VCF_filt)")
        if 'DB' not in self._tool_inputs:
            raise InvalidToolInputError("Database input is required (DB)")
        if 'BED' not in self._tool_inputs:
            raise InvalidToolInputError("AMR region input is required (BED)")
        if 'VCF_lofreq' not in self._tool_inputs:
            logger.info("Lofreq mutation input ('VCF_lofreq') was not supplied. "
                        "Low-frequency mutations will not be reported.")
        super()._check_input()

    @staticmethod
    def __get_matching_region(position: int, regions: pd.DataFrame) -> Union[dict, None]:
        """
        Returns the region matching that covers the input position (if available).
        :param position: Genome position
        :param regions: Extracted genomic regions
        :return: The matching region (if available), None otherwise
        """
        for row in regions.to_dict('records'):
            if row['start'] <= position <= row['end']:
                return row
        raise ValueError(f"Position '{position}' does not fall in AMR regions")

    @staticmethod
    def __extracts_mutation_name(row: dict, full: bool = False) -> str:
        """
        Extracts the mutation name.
        :param row: Input row
        :param full: Full name (including region name)
        :return: Mutation name
        """
        prefix = f"{row['region']['locus']}_" if full else ''
        suffix = '*' if not row['passes_filt'] else ''

        if len(row['associations']) == 0:
            effect = row['effect'] if row['effect'] is not None else 'unknown'
            return f"{prefix}{effect}{suffix}"

        unique_muts = set(x['mutation'] for x in row['associations'])
        return ';'.join(f'{prefix}{m}{suffix}' for m in unique_muts)

    def __parse_db_files(self, path_to_db: Path) -> None:
        """
        Parses the TSV files of the database.
        :param path_to_db: Path to database
        :return: None
        """
        # Parse the database with locations
        data_mut_locations = pd.read_table(path_to_db / 'mutation_locations.tsv')
        self._variant_by_key = {
            (r['position'], r['reference_nucleotide'], r['alternative_nucleotide']): r['variant']
            for r in data_mut_locations.to_dict('records')}
        logger.info(f'{len(data_mut_locations):,} mutations parsed')

        # Parse the database with AMR associations
        data_amr_association = pd.read_table(path_to_db / 'amr_associations_all.tsv')
        self._association_by_variant = collections.defaultdict(list)
        for r in data_amr_association.to_dict('records'):
            self._association_by_variant[r['variant']].append(r)
        logger.info(f'{len(data_amr_association):,} AMR associations parsed')

        # Parse antibiotics
        self._data_antibiotics = pd.read_table(path_to_db / 'antibiotics.tsv')
        self._ab_short_by_name = {r['AB']: r['Abbreviation'] for r in self._data_antibiotics.to_dict('records')}

        # Parse AMR regions
        self._data_regions = pd.read_table(self._tool_inputs['BED'][0].path, names=[
            'chr', 'start', 'end', 'locus', 'type', 'abs'])
        self._data_regions['abs_short'] = self._data_regions['abs'].apply(lambda x: ', '.join([
            self._ab_short_by_name[ab] for ab in x.split(', ')]))
        logger.info(f'{len(self._data_regions)} AMR regions parsed')

        # Parse DB version
        try:
            with open(path_to_db / 'VERSION') as handle:
                self._informs['version'] = handle.readline().strip()
        except FileNotFoundError:
            self._informs['version'] = 'n/a'

    @staticmethod
    def __parse_effect(vcf_record: Variant) -> Union[str, None]:
        """
        Parses the mutation effect from the CSQ annotation.
        Note: only extracts it for protein coding regions
        :param vcf_record: Input record
        :return: Mutation effect
        """
        # Check if BCSQ annotation is present
        bcsq = vcf_record.INFO.get('BCSQ')
        if bcsq is None:
            logger.warning(f'BCSQ info missing for: {vcf_record.CHROM}:{vcf_record.POS}')
            return None

        # Parse annotation
        parts = bcsq.split('|')
        if parts[0].startswith('&'):
            return None

        # Frameshift -> change into WHO format
        parts[0] = parts[0].lstrip('*')
        if parts[0] == 'frameshift':
            m = re.search(r'^(\d+)([A-Z])', parts[5])
            return f'p.{seq3(m.group(2))}{m.group(1)}fs' if m else 'frameshift'

        # AA change -> change into WHO format
        if parts[0] == 'missense':
            m = re.search(r'^(\d+)([A-Z])>\d+([A-Z])', parts[5])
            return f'p.{seq3(m.group(2))}{m.group(1)}{seq3(m.group(3))}' if m else 'missense'
        return parts[0]

    @staticmethod
    def __extract_af(vcf_record: Variant) -> float | None:
        """
        Extracts the allele frequency from the VCF record.
        :param vcf_record: Input VCF record
        :return: Allele frequency
        """
        af = vcf_record.INFO.get('AF')
        if af is not None:
            return af[0] if isinstance(af, (tuple, list)) else af

        # Check for DP4 (Ref-forward, Ref-reverse, Alt-forward, Alt-reverse)
        dp4 = vcf_record.INFO.get('DP4')
        if dp4 is not None:
            alt_reads = dp4[2] + dp4[3]
            total_reads = sum(dp4)
            if total_reads == 0:
                return 0.0
            return alt_reads / total_reads
        logger.info(f"Unable to extract AF from: {vcf_record}")
        return None

    def __cross_check_muts_to_db(self, vcf_input: Path, is_lofreq: bool = False) -> list[dict]:
        """
        Cross-checks the detected mutations to the database.
        :param vcf_input: Input VCF
        :param is_lofreq: is the vcf input from lofreq
        :return: List of detected AMR associations
        """
        mutations_out = []
        variants = vcfutils.parse_all_variants(vcf_input)
        logger.info(f'Parsed {len(variants):,} variants from {vcf_input}')

        for vcf_record in variants:
            effect = AMRScreen.__parse_effect(vcf_record)

            # Retrieve the corresponding region
            try:
                region = AMRScreen.__get_matching_region(vcf_record.POS, self._data_regions)
            except ValueError:
                logger.info(f'No matching region for mutation at position {vcf_record.POS}, retrying affected end')
                region = AMRScreen.__get_matching_region(vcf_record.end, self._data_regions)

            # Create the output record
            record = {
                'af': AMRScreen.__extract_af(vcf_record),
                'alt': ';'.join(str(x) for x in vcf_record.ALT),
                'effect': effect,
                'associations': [],
                'position': vcf_record.POS,
                'ref': vcf_record.REF,
                'region': region,
                'variant_type': vcf_record.var_type,
                'lofreq': is_lofreq,
                'passes_filt': False,
            }

            # Add AMR association
            for alt in vcf_record.ALT:
                # Query the database
                key = (vcf_record.POS, str(vcf_record.REF), str(alt))
                variant = self._variant_by_key.get(key)
                if variant is None:
                    continue

                # Add associations
                for association in self._association_by_variant[variant]:
                    record['associations'].append({
                        'antibiotic': association['drug'],
                        'antibiotic_short': self._ab_short_by_name[association['drug']],
                        'comment': association['comment'] if not pd.isna(association['comment']) else None,
                        'confidence': association['confidence'],
                        'effect': association['effect'],
                        'locus': association['gene'],
                        'mutation': association['mutation'],
                    })
            mutations_out.append(record)
        return mutations_out

    @staticmethod
    def __remove_bcftools_mutations_from_lofreq(input_lofreq: list[dict], input_bcftools: list[dict]) -> list[dict]:
        """
        :param input_lofreq: List of mutations identified by LoFreq
        :param input_bcftools: List of mutations identified by BCFtools
        :return: List of mutations unique to LoFreq
        """
        keys_bcf = {(m['position'], m['alt']) for m in input_bcftools}
        muts_lofreq = [m for m in input_lofreq if (m['position'], m['alt']) not in keys_bcf]
        return muts_lofreq

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Parse DB
        self.__parse_db_files(self._tool_inputs['DB'][0].path)

        # Cross-check variants in VCF with DB
        mutations_out = self.__cross_check_muts_to_db(self._tool_inputs['VCF'][0].path)
        logger.info(f"{sum(len(m['associations']) for m in mutations_out):,} AMR associations found")

        # Check if mutations passed filtering and if they were synonymous
        with VCF(str(self._tool_inputs['VCF_filt'][0].path)) as vcf_reader:
            positions_passing_filt = [v.POS for v in vcf_reader]
        mutations_out = [{
            **mut,
            'passes_filt': mut['position'] in positions_passing_filt,
        } for mut in mutations_out]
        logger.info(f"{sum(m['passes_filt'] for m in mutations_out):,}/{len(mutations_out):,} passed filtering")

        # Create text file with mutation positions
        path_bed_out = self.folder / 'mutation_positions.tsv'
        with path_bed_out.open('w') as handle:
            for pos in sorted(set(m['position'] for m in mutations_out)):
                handle.write('\t'.join(['Chromosome', str(pos)]))
                handle.write('\n')
        self._tool_outputs['TSV'] = [ToolIOFile(path_bed_out)]

        if 'VCF_lofreq' in self._tool_inputs:
            # Cross-check variants in VCF with DB separately for Lofreq variants
            mutations_lofreq = self.__cross_check_muts_to_db(self._tool_inputs['VCF_lofreq'][0].path, is_lofreq=True)
            mutations_lofreq = self.__remove_bcftools_mutations_from_lofreq(mutations_lofreq, mutations_out)
            logger.info(
                f"{sum(len(m['associations']) for m in mutations_lofreq):,} AMR associations found for Lofreq variants")

            # Include the lofreq mutations
            mutations_out.extend(mutations_lofreq)

        # Extract mutation name
        mutations_out = [{
            **row,
            'name': AMRScreen.__extracts_mutation_name(row),
            'name_full': AMRScreen.__extracts_mutation_name(row, full=True),
        } for row in mutations_out]

        # Create the JSON output file
        path_json = self.folder / 'amr_screen.json'
        with path_json.open('w') as handle:
            json.dump(mutations_out, handle, indent=2)
        logger.info(f'AMR results stored in: {path_json}')
        self._tool_outputs['JSON'] = [ToolIOFile(path_json)]
