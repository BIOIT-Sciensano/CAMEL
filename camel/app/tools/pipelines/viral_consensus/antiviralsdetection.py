import json
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
        super().__init__('antiviral detection', '0.1')

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

    def _execute_tool(self) -> None:
        """
        Executes the tool.
        :return: None
        """
        # Load the database
        subtype = self._parameters['subtype'].value
        species = self._parameters['species'].value
        data_db_muts, data_db_assoc = self.__load_db(species, subtype)

        # Parse the detected mutations
        data_detected_muts = self.__parse_detected_mutations()
        detected_muts_as_set = set(data_detected_muts[['segment', 'mutation']].itertuples(index=False, name=None))

        # Cross-check with mutations in the database
        logger.info(f'Found {len(data_db_muts):,} mutations in the database for species: {species} and subtype: {subtype}')
        detected = data_db_muts.apply(lambda x: AntiviralsDetection._is_mutation_detected(x, detected_muts_as_set), axis=1)
        logger.debug(f'Detected {sum(detected)} mutations associated with antivirals')

        associations = data_db_assoc['key'].apply(
            lambda x: AntiviralsDetection._is_association_detected(x, detected_muts_as_set, species, subtype))
        logger.debug(f'Detected {sum(associations)} associations with antivirals')

        # Create output
        path_out = self.folder / 'antivirals.json'
        with path_out.open('w') as handle:
            json.dump({
                'mutations': data_db_muts[detected].to_dict('records'),
                'associations': data_db_assoc[associations].to_dict('records')
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

    def __parse_detected_mutations(self) -> pd.DataFrame:
        """
        Parses the detected mutations.
        :return: Dataframe with detected mutations
        """
        def _parse_field(row_in: pd.Series, field: str) -> list[str]:
            value = row_in[field]
            if pd.isna(value):
                return []
            return [x.split(':')[-1] for x in value.split(',')]

        records = []
        for path_io in self._tool_inputs['TSV']:
            logger.debug(f'Parsing: {path_io.path}')
            row = pd.read_table(path_io.path).iloc[0]
            segment = path_io.path.parents[1].name.upper()

            mutations = _parse_field(row, 'aaSubstitutions')
            deletions = _parse_field(row, 'aaDeletions')
            logger.debug(
                f'Parsed {len(mutations)} mutations and {len(deletions)} deletions for {segment}')
            records += [{'mutation': m, 'segment': segment, 'type': 'AA'} for m in mutations]
            records += [{'mutation': d, 'segment': segment, 'type': 'deletion'} for d in deletions]
        return pd.DataFrame(records)
