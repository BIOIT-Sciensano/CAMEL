import pandas as pd

from pathlib import Path
from fractions import Fraction
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.core.tool import Tool


class GMats(Tool):
    """
    gMATS is the genetic Meningococcal Antigen Typing System (gMATS) for Neisseria meningitidis
    used to predict the Bexsero vaccine effectiveness.
    """

    def __init__(self) -> None:
        """
        Initializes the GMats tool.
        """
        super().__init__('GMats', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'DB' not in self._tool_inputs:
            raise InvalidToolInputError("gMATS database is required")
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError("List of input files (TSV) is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the tool.
        :return: 'None'
        """
        gmats_db = pd.read_table(self._tool_inputs['DB'][0].path, dtype=str)

        # Parse input file
        tsv_in = self._tool_inputs['TSV'][0].path
        sample_id = tsv_in.stem.replace('typing-bast-peptide-', '')
        logger.info(f'Sample ID: {sample_id}')
        antigen_info = self.__parse_input_file(tsv_in)

        # Check allele integrity
        antigen_info['Allele'] = antigen_info.apply(
            lambda x: self.__verify_integrity(x['Allele'], x['% Identity'], x['HSP/Locus length']), axis=1)

        # Determine gMATS allele
        antigen_info['Allele status'] = antigen_info.apply(
            lambda x: self.__determine_gmats_alelle(x['Locus'], x['Allele'], gmats_db), axis=1)

        # Determine global gMATS status
        status = self.__determine_gmats_status(antigen_info['Allele status'])
        antigen_info['gMATS status'] = status
        self._informs['gMATS_status'] = status

        # Get the output directory
        self.__create_output_file(antigen_info, sample_id)

    def __parse_input_file(self, file_path: Path) -> pd.DataFrame:
        """
        Parses the input file.
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

    def __determine_gmats_alelle(self, locus: str, allele: str, gmats_db: pd.DataFrame) -> str:
        """
        Determine gMATS prediction for locus/allele.
        :param locus: Locus to be assessed
        :param allele: Specific allele to be assessed
        :param gmats_db: Pandas DataFrame containing the vaccine coverage for the fHbp and NHBA loci
        :return: the gMATS prediction
        """
        fhbp_db = gmats_db.loc[gmats_db['Locus'] == 'fHbp']
        nhba_db = gmats_db.loc[gmats_db['Locus'] == 'NHBA']

        # NadA
        if locus == 'NadA_peptide':
            return 'not_covered'

        # fHbp, NHBA and PorA
        if allele == '-':
            return 'not_covered'

        if allele == '?' and locus != 'PorA_VR2':
            return 'unpredictable'

        # fHbp
        if locus == 'fHbp_peptide':
            if (fhbp_db.Allele == allele).any():
                return fhbp_db.Status.loc[fhbp_db['Allele'] == allele].iloc[0]
            else:
                return 'unpredictable'

        # NHBA
        if locus == 'NHBA_peptide':
            if (nhba_db.Allele == allele).any():
                return nhba_db.Status.loc[nhba_db['Allele'] == allele].iloc[0]
            else:
                return 'unpredictable'

        # PorA_VR2
        if locus == 'PorA_VR2':
            return 'covered' if allele == '4' else 'not_covered'

        # Other loci
        if locus == 'PorA_VR1':
            return 'not_covered'
        raise ValueError(f'Cannot determine gmats allele (locus={locus}, allele={allele})')

    def __determine_gmats_status(self, allele_status: pd.Series) -> str:
        """
        Determine global gMATS prediction for sample.
        :param allele_status: Pandas Series with the gMATS status of each allele
        :return: the global gMATS prediction
        """
        if 'covered' in list(allele_status):
            return 'covered'
        elif all(status == 'not_covered' for status in list(allele_status)):
            return 'not_covered'
        else:
            return 'unpredictable'

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
        output_file = output_dir / f'{sample_id}_gMATS.tsv'
        antigen_info.drop(['Type'], axis=1, inplace=True)
        antigen_info.to_csv(output_file, sep='\t', index=False)
        logger.info(f'Output file created: {output_file}')
        self._tool_outputs['TSV'] = [ToolIOFile(output_file)]
