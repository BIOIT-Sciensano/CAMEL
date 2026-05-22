from fractions import Fraction
from pathlib import Path

import pandas as pd
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool
from camel.app.loggers import logger


class MenDeVAR(Tool):
    """
    The MenDeVAR (Meningococcal Deduced Vaccine Antigen Reactivity) reactivity index estimates
    coverage for the Bexsero and Trumenba vaccines.
    """

    def __init__(self) -> None:
        """
        Initializes the Mendevar tool.
        """
        super().__init__('MenDeVAR', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'DB' not in self._tool_inputs:
            raise InvalidToolInputError("MenDeVAR database is required")
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError("TSV input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the tool.
        :return: None
        """
        mendevar_db = pd.read_table(self._tool_inputs['DB'][0].path, dtype=str)

        # Parse input file
        tsv_in = self._tool_inputs['TSV'][0].path
        if 'typing-bast-peptide' not in tsv_in.stem:
            logger.warning('Cannot determine sample name from TSV input, using default sample id')
            sample_id = 'sample'
        else:
            sample_id = tsv_in.stem.replace('typing-bast-peptide-', '')
        logger.info(f'Sample ID: {sample_id}')
        antigen_info = self.__parse_input_file(tsv_in)

        # Check allele integrity
        antigen_info['Allele'] = antigen_info.apply(
            lambda x: self.__verify_integrity(x['Allele'], x['% Identity'], x['HSP/Locus length']), axis=1)

        # Determine allele Bexsero index
        antigen_info['Bexsero allele status'] = antigen_info.apply(
            lambda x: self.__determine_allele_bexsero_index(x['Locus'], x['Allele'], mendevar_db), axis=1)

        # Determine Bexsero reactivity index
        bexsero_status = self.__determine_mendevar_index(antigen_info['Bexsero allele status'])
        antigen_info['Bexsero status'] = bexsero_status
        self._informs['bexsero_status'] = bexsero_status

        # Determine allele Trumenba index
        antigen_info['Trumenba allele status'] = antigen_info.apply(
            lambda x: self.__determine_allele_trumenba_index(x['Locus'], x['Allele'], mendevar_db), axis=1)

        # Determine Trumenba reactivity index
        trumenba_status = self.__determine_mendevar_index(antigen_info['Trumenba allele status'])
        antigen_info['Trumenba status'] = '-'
        antigen_info.loc[antigen_info['Locus'] == 'fHbp_peptide', 'Trumenba status'] = trumenba_status
        self._informs['trumenba_status'] = trumenba_status

        # Get the output directory
        self.__create_output_file(antigen_info, sample_id)

    def __parse_input_file(self, file_path: Path) -> pd.DataFrame:
        """
        Parses the input file for MenDeVAR coverage prediction.
        :param file_path: Path to input text file
        :return: DataFrame with allele info
        """
        input_file = pd.read_table(file_path, dtype=str)
        input_file = input_file[input_file['Locus'] != 'PorA_VR1']
        return input_file

    def __verify_integrity(self, allele: str, identity: str, coverage: str) -> str:
        """
        Verifies that the allele is a perfect hit i.e. that it shows 100% identity and 100% coverage.
        :param allele: Allele
        :param identity: Allele identity %
        :param coverage: Allele coverage length
        :return: Allele or allele marked with an asterisk(*) in case of imperfect match
        """
        # Skip undetermined alleles
        if allele in ['-', '?']:
            return allele

        # Verify assigned alleles
        if identity == '100.00' and int(Fraction(coverage)) == 1:
            return allele
        else:
            return allele + '*'

    def __determine_allele_bexsero_index(self, locus: str, allele: str, mendevar_db: pd.DataFrame) -> str:
        """
        Determine Bexsero reactivity index for locus/allele.
        :param locus: Locus to be assessed
        :param allele: Specific allele to be assessed
        :param mendevar_db: Pandas DataFrame containing the vaccine coverage for the fHbp, NHBA, NadA and PorA_VR2 loci
        :return: Bexsero reactivity indexes
        """
        fhbp_db = mendevar_db.loc[mendevar_db['Locus'] == 'fHbp']
        nhba_db = mendevar_db.loc[mendevar_db['Locus'] == 'NHBA']
        nada_db = mendevar_db.loc[mendevar_db['Locus'] == 'NadA']

        # Missing allele
        if allele == '-':
            return 'none'

        # New allele
        if allele == '?':
            return 'insufficient data'

        # fHbp
        if locus == 'fHbp_peptide':
            if (fhbp_db.Allele == allele).any():
                return fhbp_db.Bexsero_status.loc[fhbp_db['Allele'] == allele].iloc[0]
            else:
                return 'insufficient data'

        # NHBA
        if locus == 'NHBA_peptide':
            if (nhba_db.Allele == allele).any():
                return nhba_db.Bexsero_status.loc[nhba_db['Allele'] == allele].iloc[0]
            else:
                return 'insufficient data'

        # NadA
        if locus == 'NadA_peptide':
            if (nada_db.Allele == allele).any():
                return nada_db.Bexsero_status.loc[nada_db['Allele'] == allele].iloc[0]
            else:
                return 'insufficient data'

        # PorA_VR2
        if locus == 'PorA_VR2':
            return 'exact match' if allele == '4' else 'none'

        # Other loci
        if locus == 'PorA_VR1':
            return 'none'
        raise ValueError(f'Cannot determine allele Bexsero reactivity index (locus={locus}, allele={allele})')

    def __determine_allele_trumenba_index(self, locus: str, allele: str, mendevar_db: pd.DataFrame) -> str:
        """
        Determine Trumenba reactivity index for locus/allele.
        :param locus: Locus to be assessed
        :param allele: Specific allele to be assessed
        :param mendevar_db: Pandas DataFrame containing the vaccine coverage for the fHbp locus
        :return: Trumenba reactivity indexes
        """
        fhbp_db = mendevar_db.loc[mendevar_db['Locus'] == 'fHbp']

        # fHbp
        if locus == 'fHbp_peptide':
            # Missing allele
            if allele == '-':
                return 'none'

            # New allele
            if allele == '?':
                return 'insufficient data'

            # Known allele
            if (fhbp_db.Allele == allele).any():
                return fhbp_db.Trumenba_status.loc[fhbp_db['Allele'] == allele].iloc[0]
            else:
                return 'insufficient data'

        # Other loci
        if locus in ('NHBA_peptide', 'NadA_peptide', 'PorA_VR2', 'PorA_VR1'):
            return 'NA'
        raise ValueError(f'Cannot determine allele Trumenba reactivity index (locus={locus}, allele={allele})')

    def __determine_mendevar_index(self, allele_status: pd.Series) -> str:
        """
        Determine global MenDeVAR reactivity index for the sample.
        :param allele_status: Pandas Series with the reactivity index of each allele
        :return: the global MenDeVAR reactivity index
        """
        if 'exact match' in list(allele_status):
            return 'exact match'
        elif 'cross_reactive' in list(allele_status):
            return 'cross_reactive'
        elif 'insufficient data' in list(allele_status):
            return 'insufficient data'
        else:
            return 'none'

    def __create_output_file(self, antigen_info: pd.DataFrame, sample_id: str) -> None:
        """
        Write output tsv.
        :param antigen_info: Final antigen_info Pandas DataFrame
        :return: None
        """
        # Get the output directory
        if 'output_directory' in self._parameters:
            output_dir = Path(self._parameters["output_directory"].value)
        else:
            output_dir = self.folder

        # Create output file
        output_file = output_dir / f'{sample_id}_mendevar.tsv'
        antigen_info.drop(['Type'], axis=1, inplace=True)
        antigen_info.to_csv(output_file, sep='\t', index=False)
        logger.info(f'Output file created: {output_file}')
        self._tool_outputs['TSV'] = [ToolIOFile(output_file)]
