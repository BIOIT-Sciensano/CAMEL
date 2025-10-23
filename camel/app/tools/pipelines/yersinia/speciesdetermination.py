import pandas as pd
from pathlib import Path

from camel.app.core.tool import Tool
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.errors import InvalidToolInputError


class SpeciesDetermination(Tool):
    """
    This tool is used to determine the species and lineage of Yersinia isolates.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Yersinia: species determination', '1.0')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Find species matching to sequence types and over the species-specific threshold
        taxonomic = pd.read_table(self._tool_inputs['TSV_taxonomic'][0].path)
        profile_matches = pd.read_table(self._tool_inputs['TSV_profile_matches'][0].path, dtype={'ST':int, 'proportion_match':float})
        output = self.__find_matches(taxonomic, profile_matches)

        # Save output data
        output.to_csv(self._folder / self._parameters['output_filename'].value, sep="\t", index=False)
        self._tool_outputs['TSV'] = [ToolIOFile(self.folder / Path(self._parameters['output_filename'].value))]
        if output.shape[0] == 0:
            self._informs.update({'found_match': False})
        else:
            best_match = output.loc[0]
            self._informs.update({'found_match': True, 'best_match': dict(best_match)})

    def _check_input(self) -> None:
        """
        Check input for species determination tool.
        :return: None
        """
        if 'TSV_taxonomic' not in self._tool_inputs:
            raise InvalidToolInputError("TSV taxonomic file is required")
        if 'TSV_profile_matches' not in self._tool_inputs:
            raise InvalidToolInputError("TSV profile matches file is required")
        super()._check_input()

    def __find_matches(self, taxonomic: pd.DataFrame, profile_matches: pd.DataFrame) -> pd.DataFrame:
        """
        Finds sequence type matches that are over the thresholds set for the species and lineage.
        :param taxonomic: dataframe with correspondence ST, species, and thresholds
        :param profile_matches: dataframe with ST profiles and proportion of matching loci found during typing
        :return: dataframe with species matching over their threshold
        """
        result = pd.DataFrame(columns=['species', 'lineage', 'biotype', 'serotype', 'match', 'threshold'])
        for _, row in taxonomic.iterrows():
            # Find proportion of matching loci in typing result.
            match = profile_matches.loc[profile_matches['ST'] == row['cgST'], 'proportion_match'].values[0]
            if match >= row['threshold']:
                # Check if the lineage is not already in the results
                i = result[(result['species'] == row['species']) & (result['lineage'] == row['lineage'])].index
                if len(i) > 0:
                    # If the lineage is already in the results, only report the match with the highest proportion
                    # of matching loci
                    result.loc[i[0], 'match'] = max(match, result.loc[i[0], 'match'])
                else:
                    # If the lineage is not in the results, add a row
                    result.loc[len(result.index)] = [row['species'], row['lineage'], row['biotype'],
                                                     row['serotype'], match, row['threshold']]
        result = result.fillna('NA')
        return result.sort_values(by='match', ascending=False, ignore_index=True)
