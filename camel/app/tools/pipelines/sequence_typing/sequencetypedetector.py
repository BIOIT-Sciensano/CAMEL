from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

from camel.app.toolkits.sequencetyping.typinghitbase import TypingHitBase
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.core.tool import Tool


@dataclass(frozen=True, unsafe_hash=True)
class STProfile:
    """
    Dataclass to hold ST profile data.
    """
    name: str
    metadata: list[tuple[str, str]] = field(hash=False)
    alleles: Optional[dict[str, str]] = field(hash=False)


class SequenceTypeDetector(Tool):
    """
    Tool that manages MLST schemes. Also reports scheme metadata information in the informs.
    """

    SYMBOL_NO_ST = 'ND'

    def __init__(self):
        """
        Initialize this tool.
                :return: None
        """
        super().__init__('Typing: Sequence Type Detector', '0.1')
        self._wildcards = None
        self._symbol_allele_absent = None
        self._metadata_columns = []

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Parse input data
        hit_by_locus = {hit.locus: hit for hit in [
            io.value for io in self._tool_inputs['hits_nucl'] + self._tool_inputs['hits_pept']]}
        profiles = self.__parse_profiles(list(hit_by_locus.keys()))

        # Determine best profile
        nb_matches_by_profile = self.__get_nb_matches_by_profile(profiles, hit_by_locus)
        best_profile = sorted(nb_matches_by_profile.items(), key=lambda x: -x[-1])[0][0]
        percent_matching = 100 * nb_matches_by_profile[best_profile] / len(best_profile.alleles)
        is_detected = percent_matching >= int(self._parameters['min_percent_detected'].value)

        # Write all profiles to a tsv file
        if ('write_tsv' in self._parameters.keys()) and self._parameters['write_tsv'].value:
            num_alleles = len(list(nb_matches_by_profile.keys())[0].alleles)
            all_matches = pd.DataFrame([(
                profile.name, hits / num_alleles) for profile, hits in nb_matches_by_profile.items()],
                columns=['ST', 'proportion_match'])
            all_matches.to_csv(self._folder / Path(self._parameters['output_filename'].value), sep="\t", index=False)
            self._tool_outputs['TSV'] = [ToolIOFile(self._folder / Path(self._parameters['output_filename'].value))]

        # Save output data
        self._informs.update({
            'is_detected': is_detected,
            'nb_detected': nb_matches_by_profile[best_profile],
            'nb_loci': len(best_profile.alleles),
            'percent_detected': percent_matching,
            'symbol': best_profile.name if is_detected else SequenceTypeDetector.SYMBOL_NO_ST,
            'metadata': [(k, v if is_detected else '-') for k, v in best_profile.metadata]
        })

    def _check_input(self) -> None:
        """
        Checks whether the input is correct.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError("Sequence type profiles ('TSV') are required.")
        if len(self._tool_inputs['hits_nucl']) + len(self._tool_inputs['hits_pept']) == 0:
            raise InvalidToolInputError("Typing hits are required.")
        super()._check_input()

    def __parse_profiles(self, gene_names: list[str]) -> list[STProfile]:
        """
        Parses the sequence type profiles.
        :param gene_names: Name of the genes
        :return: List of profiles
        """
        # Parse input data
        data_in = pd.read_table(self._tool_inputs['TSV'][0].path, dtype=str)
        cols_metadata = [c for c in data_in.columns if c not in gene_names]
        logger.info(f'Metadata columns: {cols_metadata}')
        cols_alleles = [c for c in data_in.columns if c in gene_names]
        logger.info(f'Gene columns: {cols_alleles}')

        # Construct the profiles
        profiles = []
        for row in data_in.fillna('n/a').to_dict('records'):
            profiles.append(STProfile(
                name=row[str(data_in.columns[0])],
                alleles={c: row[c] for c in cols_alleles},
                metadata=[(c, row[c] if not pd.isna(row[c]) else '-') for c in cols_metadata],
            ))
        logger.info(f'Parsed {len(profiles):,} profiles')
        return profiles

    def __alleles_match(self, detected_hit: TypingHitBase, profile_allele: str) -> bool:
        """
        Checks whether two alleles match.
        :param detected_hit: Detected allele
        :param profile_allele: Allele from the ST profile
        :return: True if the alleles match
        """
        if profile_allele == self._symbol_allele_absent:
            return detected_hit.allele_id == '-'
        elif profile_allele in self._wildcards:
            return True
        return (detected_hit.is_perfect_hit()) and (detected_hit.allele_id == profile_allele)

    def __get_nb_matches_by_profile(
            self, profiles: list[STProfile], hit_by_locus: dict[str, TypingHitBase]) -> dict:
        """
        Returns the number of matches per profile.
        :param profiles: Sequence type profiles
        :param hit_by_locus: Alleles for the genes.
        :return: Nb. of matches by profile
        """
        if 'allele_absent_symbol' in self._parameters:
            self._symbol_allele_absent = self._parameters['allele_absent_symbol'].value
        self._wildcards = self._parameters['allele_wildcards'].value.split(',')
        logger.debug("Wildcards: [{}]¸ Symbol allele absent: {}".format(
            ', '.join(self._wildcards), self._symbol_allele_absent))

        nb_matches_by_profile = {}
        for profile in profiles:
            nb_matching = sum([
                self.__alleles_match(hit, profile.alleles[gene_name]) for gene_name, hit in hit_by_locus.items()])
            nb_matches_by_profile[profile] = nb_matching
        return nb_matches_by_profile
