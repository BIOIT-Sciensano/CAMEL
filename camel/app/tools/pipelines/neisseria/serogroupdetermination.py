
from enum import Enum

from camel.app.toolkits.sequencetyping.typingblasthit import TypingBlastHit
from camel.app.core.tool import Tool


class DetectionCategory(Enum):
    ALL_PERFECT = 1
    IMPERFECT_IDENTITY = 2
    IMPERFECT_ALL = 3
    OTHER = 4


class SerogroupDetermination(Tool):
    """
    This tool is used to determine the serogroup of Neisseria meningitidis isolates.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Neisseria: serogroup determination', '0.1')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        serogroup_stats = []

        for key in self._tool_inputs:
            # Extract hits and sort them
            hits = [h.value for h in self._tool_inputs[key]]
            hits.sort(key=lambda h: h.locus)

            # Determine the number of perfect hits
            nb_detected = len([h for h in hits if h.allele_id != '-'])
            nb_perfect_hits = len([h for h in hits if h.is_perfect_hit()])
            serogroup_stats.append({
                'name': key.split('_')[-1].upper(),
                'category': self.__get_detection_category(hits).value,
                'nb_loci_total': len(hits),
                'nb_hits': nb_detected,
                'nb_hits_perfect': nb_perfect_hits,
                'fraction_detected': nb_detected / len(hits),
                'fraction_detected_perfect': nb_perfect_hits / len(hits),
                'color_per_hit': [(h.locus, h.color) for h in hits]
            })

        # Sort serogroups based on category and fraction of perfect hits
        serogroup_stats.sort(key=lambda x: (x['category'], -x['fraction_detected'], -x['fraction_detected_perfect']))
        self._informs['serogroups_sorted'] = serogroup_stats
        self._informs['detected_serogroup'] = serogroup_stats[0]['name'] if \
            serogroup_stats[0]['fraction_detected'] >= 0.6 else 'NA'

    def __get_detection_category(self, hits: list[TypingBlastHit]) -> DetectionCategory:
        """
        Returns the category of hits for the given list of hits.
        :param hits: Hits
        :return: Detection type
        """
        if all(h.is_perfect_hit() for h in hits):
            return DetectionCategory.ALL_PERFECT
        elif all(h.is_full_length() for h in hits):
            return DetectionCategory.IMPERFECT_IDENTITY
        elif not any(h.allele_id == '-' for h in hits):
            return DetectionCategory.IMPERFECT_ALL
        return DetectionCategory.OTHER
