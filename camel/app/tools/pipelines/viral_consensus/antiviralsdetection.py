import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool
from camel.app.loggers import logger


class AntiviralsDetection(Tool):
    """
    Detects antiviral mutations by cross-checking the Nextclade output with an in-house database.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('antiviral detection', '0.2')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'DB' not in self._tool_inputs:
            raise InvalidToolInputError('DB input is required')
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError('TSV input from Nextclade is required')
        super()._check_input()

    @staticmethod
    def _is_mutation_detected(db_row: pd.Series, detected_muts: set[tuple[str, str]]) -> bool:
        """
        Checks if the target mutation from the DB is detected in this sample.
        :param db_row: Database row
        :param detected_muts: List of detected mutations
        :return: True if detected, False otherwise
        """
        return (db_row['segment'], db_row['mutation']) in detected_muts

    @staticmethod
    def _is_association_detected(mut_string: str, detected: set[tuple[str, str]], species: str, subtype: str) -> bool:
        """
        Returns True if the target mutation was detected, False otherwise.
        :param mut_string: Target mutation string
        :param detected: List of detected mutations
        :param species: Species name
        :param subtype: Subtype name
        """
        for m in mut_string.split('+'):
            species_, subtype_, segment, mut = m.split('_')
            if species_ != species or subtype_ != subtype:
                return False
            if (segment, mut) not in detected:
                return False
        return True

    @staticmethod
    def parse_position(mutation: str) -> int | None:
        """
        Extracts the amino acid position from a single mutation (e.g. S211N or N245-).
        :param mutation: Mutation string
        :return: Position, None if the mutation does not target a single position (e.g. combinations)
        """
        m = re.fullmatch(r'[A-Z*](\d+)[A-Z*-]', mutation)
        return int(m.group(1)) if m else None

    @staticmethod
    def _is_mutation_assessable(
            db_row: pd.Series, cds_covered: set[str], unknown_ranges: dict[str, list[tuple[int, int]]]) -> bool:
        """
        Checks if the position of the target mutation from the DB is covered by the consensus sequence.
        :param db_row: Database row
        :param cds_covered: CDSs covered by the consensus sequence
        :param unknown_ranges: Unknown amino acid ranges by CDS
        :return: True if the mutation can be assessed, False otherwise
        """
        position = AntiviralsDetection.parse_position(db_row['mutation'])
        if position is None:
            # Combinations are assessed through their individual mutations
            return True
        if db_row['segment'] not in cds_covered:
            return False
        return not any(start <= position <= end for start, end in unknown_ranges.get(db_row['segment'], []))

    @staticmethod
    def _split_cds_values(value: str, cds_default: str) -> list[tuple[str, str]]:
        """
        Splits a comma-separated Nextclade field with CDS prefixes (e.g. 'F:S211N,G:A74V').
        :param value: Field value
        :param cds_default: CDS name used when a value has no CDS prefix
        :return: List of (CDS, value) tuples
        """
        if not value:
            return []
        records = []
        for item in value.split(','):
            cds, _, val = item.rpartition(':')
            records.append(((cds or cds_default).upper(), val))
        return records

    def _execute_tool(self) -> None:
        """
        Executes the tool.
        :return: None
        """
        # Load the database
        subtype = self.get_param_value('subtype')
        species = self.get_param_value('species')
        data_db_muts, data_db_assoc = self.__load_db(species, subtype)

        # Parse the Nextclade output
        detected_muts, cds_covered, unknown_ranges = self.__parse_nextclade_output()

        # Cross-check with mutations in the database
        logger.info(f'Found {len(data_db_muts):,} mutations in the database for species: {species} and subtype: {subtype}')
        detected = pd.Series(
            [AntiviralsDetection._is_mutation_detected(row, detected_muts) for _, row in data_db_muts.iterrows()],
            index=data_db_muts.index, dtype=bool)
        logger.debug(f'Detected {sum(detected)} mutations associated with antivirals')

        associations = pd.Series(
            [AntiviralsDetection._is_association_detected(key, detected_muts, species, subtype)
             for key in data_db_assoc['key']],
            index=data_db_assoc.index, dtype=bool)
        logger.debug(f'Detected {sum(associations)} associations with antivirals')

        # Check which database mutations cannot be assessed (e.g. low coverage)
        assessable = pd.Series(
            [AntiviralsDetection._is_mutation_assessable(row, cds_covered, unknown_ranges)
             for _, row in data_db_muts.iterrows()],
            index=data_db_muts.index, dtype=bool)
        if (~assessable).any():
            logger.warning(f'{sum(~assessable)} database mutations cannot be assessed due to missing sequence data')

        # Create output
        path_out = self.folder / 'antivirals.json'
        with path_out.open('w') as handle:
            json.dump({
                'mutations': data_db_muts[detected].to_dict('records'),
                'associations': data_db_assoc[associations].to_dict('records'),
                'mutations_not_assessable': data_db_muts[~assessable].to_dict('records'),
            }, handle, indent=2)
        self._tool_outputs['JSON'] = [ToolIOFile(path_out)]

    def __load_db(self, species: str, subtype: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Loads the database and filters by species and subtype.
        :param species: Detected species
        :param subtype: Detected subtype
        :return: Database mutations, database associations
        """
        # Mutations
        path_tsv = Path(self._tool_inputs['DB'][0].path, 'mutations.tsv')
        if not path_tsv.exists():
            raise FileNotFoundError(f'Mutations file not found: {path_tsv}')
        data_db_mutations = pd.read_table(path_tsv, keep_default_na=False, na_values='-')
        data_db_mutations = data_db_mutations[
            (data_db_mutations['species'] == species) & (data_db_mutations['subtype'] == subtype)]
        if len(data_db_mutations) == 0:
            logger.warning(f'No mutations in database for subtype: {subtype}')

        # Associations
        path_tsv_assoc = Path(self._tool_inputs['DB'][0].path, 'associations.tsv')
        if not path_tsv_assoc.exists():
            raise FileNotFoundError(f'Associations file not found: {path_tsv_assoc}')
        data_db_associations = pd.read_table(path_tsv_assoc, keep_default_na=False, na_values='-')
        data_db_associations = data_db_associations[
            data_db_associations['key'].str.startswith(f'{species}_{subtype}_')]
        return data_db_mutations, data_db_associations

    def __parse_nextclade_output(self) -> tuple[set[tuple[str, str]], set[str], dict[str, list[tuple[int, int]]]]:
        """
        Parses the detected mutations and the coverage information from the Nextclade output.
        :return: Detected mutations as (CDS, mutation), covered CDSs, unknown amino acid ranges by CDS
        """
        detected_muts = set()
        cds_covered = set()
        unknown_ranges = defaultdict(list)
        for path_io in self._tool_inputs['TSV']:
            logger.debug(f'Parsing: {path_io.path}')
            # Default NA values are disabled, as 'NA' is a valid CDS name
            row = pd.read_table(path_io.path, keep_default_na=False, dtype=str).iloc[0]
            segment = path_io.path.parents[1].name.upper()

            # Mutations
            mutations = AntiviralsDetection._split_cds_values(row.get('aaSubstitutions', ''), segment)
            deletions = AntiviralsDetection._split_cds_values(row.get('aaDeletions', ''), segment)
            logger.debug(f'Parsed {len(mutations)} mutations and {len(deletions)} deletions for {segment}')
            detected_muts.update(mutations + deletions)

            # Coverage
            cds_failed = {cds.strip().upper() for cds in row.get('failedCdses', '').split(',') if cds.strip()}
            for cds, coverage in AntiviralsDetection._split_cds_values(row.get('cdsCoverage', ''), segment):
                if (float(coverage) > 0) and (cds not in cds_failed):
                    cds_covered.add(cds)
            for cds, range_ in AntiviralsDetection._split_cds_values(row.get('unknownAaRanges', ''), segment):
                start, _, end = range_.partition('-')
                unknown_ranges[cds].append((int(start), int(end or start)))
        return detected_muts, cds_covered, dict(unknown_ranges)
