import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

from camelcore.app.io.tooliovalue import ToolIOValue

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool
from camel.app.loggers import logger
from camel.app.toolkits.blast.blastformat7parser import BlastFormat7Parser


@dataclass
class SpaTypingHit:
    """
    This class represents a hit for a spa type.
    """
    spa_type: str
    repeats: list[int]
    length: int
    percent_identity: float
    percent_covered: float
    blast_output: dict[str, Any]

    def is_perfect(self) -> bool:
        """
        Returns true if it is a perfect hit.
        :return: True if perfect hit
        """
        return self.percent_identity == 100.0 and self.percent_covered == 100.0

    @property
    def genomic_coordinates(self) -> str:
        """
        Returns the genomic coordinates for the hit.
        :return: Genomic coordinates
        """
        return f"{self.blast_output['qseqid']}:{self.blast_output['qstart']}-{self.blast_output['qend']}"

    @property
    def strand(self) -> str:
        """
        Returns the strand on which the spa sequence is located.
        :return: Strand
        """
        return self.blast_output['sstrand']


class SpaTyping(Tool):
    """
    This class is used to perform spa-typing based on blastn output.
    """

    BLASTN_OUTPUT_FORMAT = '"7 qseqid sseqid pident slen qstart qend length sseq sstrand"'

    def __init__(self):
        """
        Initializes this tool.
        """
        super().__init__('Spa typing', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError('Tabular BLAST input is required')
        if 'CSV_profiles' not in self._tool_inputs:
            raise InvalidToolInputError('spa type profiles input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        profiles = SpaTyping.__parse_spa_type_profiles(self._tool_inputs['CSV_profiles'][0].path)
        filtered_blast_output = self.__parse_blast_output(self._tool_inputs['TSV'][0].path)
        spa_type_hits = self.__convert_blast_output(filtered_blast_output, profiles)
        self.__set_output(spa_type_hits, profiles)

    @staticmethod
    def __parse_spa_type_profiles(profiles_tsv: Path) -> dict[str, list[int]]:
        """
        Parses the spa types from the tabular profiles file.
        :return: Dictionary of spa types.
        """
        profiles = {}
        with open(profiles_tsv) as handle:
            for line in handle.readlines():
                if line.startswith('NT'):
                    continue
                profile, repeats = line.strip().split(',')
                profiles[profile] = [int(x) for x in repeats.split('-')]
        return profiles

    def __parse_blast_output(self, blast_output_path: Path) -> list[dict]:
        """
        Parses the tabular blast output, returns a list of hits sorted by percent covered and percent identity.
        Hits that cover less than 90% or have less than 90% identity are filtered out.
        :param blast_output_path: Path to the BLAST output file
        :return: List of hits
        """
        blast_out = BlastFormat7Parser.parse_output_file(blast_output_path)
        filtered_hits = []
        for h in blast_out:
            if (not SpaTyping.__calculate_percent_covered(h) > 90) or (not h['pident'] > 90):
                continue
            filtered_hits.append(h)
        filtered_hits.sort(key=lambda x: self.__calculate_percent_covered(x), reverse=True)
        return filtered_hits

    def __convert_blast_output(self, blast_output: list[dict[str, Any]], profiles: dict[str, list[int]]) -> \
            list[SpaTypingHit]:
        """
        Converts the BLAST output to spa typing hits.
        :param blast_output: BLAST output
        :return: List of spa typing hits
        """
        hits = []
        for h in blast_output[:10]:
            percent_covered = self.__calculate_percent_covered(h)
            hits.append(SpaTypingHit(h['sseqid'], profiles[h['sseqid']], h['slen'], h['pident'], percent_covered, h))
        hits.sort(reverse=True, key=lambda x: (x.percent_covered, x.percent_identity))
        logger.debug(f"{len(hits)} hits detected")
        return hits

    def __set_output(self, spa_type_hits: list[SpaTypingHit], profiles: dict[str, list[int]]) -> None:
        """
        Sets the output of this tool.
        :param spa_type_hits: List of spa type hits
        :param profiles: Dictionary of profiles with associated repeats
        :return: None
        """
        perfect_hits = [h for h in spa_type_hits if h.is_perfect()]
        logger.debug(f"{len(perfect_hits)} perfect hits found")
        if len(perfect_hits) == 1:
            self._informs['spa_type'] = perfect_hits[0].spa_type
            self._informs['genomic_coordinates'] = perfect_hits[0].genomic_coordinates
            self._informs['strand'] = perfect_hits[0].strand
        elif len(perfect_hits) > 1:
            all_types = [x.spa_type for x in perfect_hits]
            self._informs['spa_type'] = f"ambiguous ({', '.join(sorted(all_types))})"
        else:
            self._informs['spa_type'] = 'NA'
        self._informs['spa_type_repeats'] = profiles.get(self._informs['spa_type'])
        self._informs['db_info'] = SpaTyping.__get_db_info(self._tool_inputs['CSV_profiles'][0].path.parent)
        self._tool_outputs['VAL_hits'] = [ToolIOValue(h) for h in spa_type_hits]

    @staticmethod
    def __calculate_percent_covered(blast_data: dict[str, Any]) -> float:
        """
        Calculates the percentage of the subject sequence that is covered by the alignment.
        :param blast_data: Dictionary containing the BLAST output data
        :return: Percentage covered
        """
        return min(100 * blast_data['length'] / blast_data['slen'], 100.0)

    @staticmethod
    def __get_db_info(dir_db: Path) -> Optional[dict]:
        """
        Returns the date of the last database update.
        :param dir_db: Database directory
        :return: Database information
        """
        path_json = dir_db / 'db_update_info.json'
        if not path_json.exists():
            logger.info(f'No database update file found: {path_json}')
            return None
        with path_json.open() as handle:
            return json.load(handle)
