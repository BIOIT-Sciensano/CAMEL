import json
from pathlib import Path
from typing import Union

import pandas as pd
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.tool import Tool
from camel.app.loggers import logger
from camel.app.toolkits.mycobacterium import amrutils
from camel.app.toolkits.mycobacterium.amrutils import ConfidenceLevel


class AMRPhenotypePredictor(Tool):
    """
    This tool determines the AMR type of based on the detected mutations.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Mycobacterium: AMR phenotype predictor', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'DIR_DB' not in self._tool_inputs:
            raise InvalidToolInputError("Database directory is required (DIR_DB).")
        if 'JSON' not in self._tool_inputs:
            raise InvalidToolInputError("AMR association information is required (JSON).")
        super()._check_input()

    @staticmethod
    def __get_mutations_by_antibiotic(input_file: Path) -> dict[str, dict[str, dict]]:
        """
        Returns the detected mutations grouped by antibiotic and confidence level.
        :param input_file: Input file with mutation information
        :return: Mutations grouped by antibiotic
        """
        mutations_by_ab = {}
        with input_file.open() as handle:
            data_mutations = json.load(handle)

        for mutation in data_mutations:
            if mutation['lofreq'] is True:
                continue
            for association in mutation['associations']:
                if association['antibiotic'] not in mutations_by_ab:
                    mutations_by_ab[association['antibiotic']] = {}
                confidence = ConfidenceLevel(association['confidence'])
                if confidence.value not in mutations_by_ab[association['antibiotic']]:
                    mutations_by_ab[association['antibiotic']][confidence.value] = []
                mutations_by_ab[association['antibiotic']][confidence.value].append(mutation)
        return mutations_by_ab

    @staticmethod
    def __predict_phenotype(mutations_by_confidence: Union[dict[amrutils.ConfidenceLevel, list], None]) -> str:
        """
        Predicts the phenotype based on the detected mutations.
        :param mutations_by_confidence: Mutations group by confidence level.
        :return: Predicted phenotype
        """
        if mutations_by_confidence is None:
            return 'S'
        confidence_levels = [amrutils.ConfidenceLevel(x) for x in mutations_by_confidence.keys()]
        if amrutils.ConfidenceLevel.ASSOC_R in confidence_levels:
            return 'R'
        elif amrutils.ConfidenceLevel.ASSOC_R_int in confidence_levels:
            return 'R (int.)'
        elif amrutils.ConfidenceLevel.ASSOC_S_int in confidence_levels:
            return 'S (int.)'
        return 'S'

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Parse input data
        data_ab = pd.read_table(Path(self._tool_inputs['DIR_DB'][0].path, 'antibiotics.tsv'))
        mutations_by_ab = self.__get_mutations_by_antibiotic(Path(self._tool_inputs['JSON'][0].path))

        # Create output data structure
        data_out = []
        for data_ab in data_ab.to_dict('records'):
            data_out.append({
                'category': data_ab['Category'],
                'name': data_ab['AB'],
                'abbreviation': data_ab['Abbreviation'],
                'mutations': mutations_by_ab.pop(data_ab['AB'], {})
            })

        if len(mutations_by_ab) > 0:
            logger.error(f'Unparsed mutations for: {mutations_by_ab.keys()}')
            raise ToolExecutionError(self.name, 'unparsed mutations')

        # Predict phenotype
        for row in data_out:
            row['phenotype'] = AMRPhenotypePredictor.__predict_phenotype(row['mutations'])

        # Save output file
        output_path = self.folder / 'mutations_by_ab.json'
        with output_path.open('w') as handle:
            json.dump(data_out, handle, indent=2)
        self._tool_outputs['JSON'] = [ToolIOFile(output_path)]
