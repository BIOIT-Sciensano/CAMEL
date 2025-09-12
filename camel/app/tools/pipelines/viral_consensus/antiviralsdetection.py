import json
from pathlib import Path

import pandas as pd

from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


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

    def _execute_tool(self) -> None:
        """
        Executes the tool.
        :return: None
        """
        # Load the database
        subtype = self._parameters['subtype'].value
        species = self._parameters['species'].value
        data_db_muts, data_db_assoc = self.__load_db(subtype)

        # Parse the detected mutations
        data_detected_muts = self.__parse_detected_mutations()
        if not data_detected_muts.empty:
            detected_muts_as_tuple = list(data_detected_muts[['segment', 'mutation']].itertuples(index=False, name=None))
        else:
            detected_muts_as_tuple = []

        # Cross-check with mutations in the database
        detected = data_db_muts.apply(lambda x: (x['segment'], x['mutation']) in detected_muts_as_tuple, axis=1)
        logger.debug(f'Detected {sum(detected)} mutations associated with antivirals')
        associations = data_db_assoc['key'].apply(
            lambda x: AntiviralsDetection.is_detected(x, detected_muts_as_tuple, species, subtype))
        logger.debug(f'Detected {sum(associations)} associations with antivirals')

        # Create output
        path_out = self.folder / 'antivirals.json'
        with path_out.open('w') as handle:
            json.dump({
                'mutations': data_db_muts[detected].to_dict('records'),
                'associations': data_db_assoc[associations].to_dict('records')
            }, handle, indent=2)
        self._tool_outputs['JSON'] = [ToolIOFile(path_out)]

    def __load_db(self, subtype: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Loads the database.
        :param subtype: Detected subtype
        :return: Database mutations, database associations
        """
        # Mutations
        path_tsv = Path(self._tool_inputs['DB'][0].path, 'mutations.tsv')
        if not path_tsv.exists():
            raise FileNotFoundError(f'Mutations file not found: {path_tsv}')
        data_db_mutations = pd.read_table(path_tsv, keep_default_na=False, na_values='-')
        data_db_mutations = data_db_mutations[data_db_mutations['subtype'] == subtype]
        if len(data_db_mutations) == 0:
            logger.warning(f'No mutations in database for subtype: {len(subtype)}')
        logger.info(f'Mutations for subtype: {len(data_db_mutations):,}')

        # Associations
        path_tsv_assoc = Path(self._tool_inputs['DB'][0].path, 'associations.tsv')
        if not path_tsv_assoc.exists():
            raise FileNotFoundError(f'Associations file not found: {path_tsv}')
        data_db_associations = pd.read_table(path_tsv_assoc, keep_default_na=False, na_values='-')
        return data_db_mutations, data_db_associations

    def __parse_detected_mutations(self) -> pd.DataFrame:
        """
        Parses the detected mutations.
        :return: Dataframe with detected mutations
        """
        mutations_all = []
        for path_io in self._tool_inputs['TSV']:
            logger.debug(f'Parsing: {path_io.path}')
            data_n3 = pd.read_table(path_io.path)
            if not pd.isna(data_n3.loc[0, 'aaSubstitutions']):
                mutations_segment = data_n3.loc[0, 'aaSubstitutions'].split(',')
            else:
                mutations_segment = []

            # Remove leading segment names
            mutations_segment = [x.split(':')[-1] for x in mutations_segment]
            logger.debug(f'Parsed: {len(mutations_segment)} mutations from: {path_io.path.name}')
            data_segment = pd.DataFrame({
                'mutation': mutations_segment,
                'segment': path_io.path.parent.name.upper()
            })
            mutations_all.extend(data_segment.to_dict('records'))
        return pd.DataFrame(mutations_all)

    @staticmethod
    def is_detected(mut_string: str, detected: list[tuple[str, str]], species: str, subtype: str) -> bool:
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
