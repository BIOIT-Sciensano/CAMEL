import pandas as pd

from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool
from camel.app.core.utils import toolutils
from camel.app.toolkits.genedetection.genedetectionmixedhit import GeneDetectionMixedHit


class MixedHitAggregatingAndFiltering(Tool):
    """
    Performs mixed hit filtering (KMA + BLAST).
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('Gene Detection: Mixed Hit Filtering', '0.1')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        hits_combined = self.__combine_hits()
        summary_tables = self.__create_summary_tables(hits_combined)
        self._tool_outputs['VAL_hits'] = [ToolIOValue(hit) for hit in hits_combined]
        self._tool_outputs['TSV_summary'] = [ToolIOValue(x) for x in summary_tables]

    def _check_input(self) -> None:
        """
        Checks whether the input is correct.
        :return: None
        """
        toolutils.check_input(self, keys_required=['KMA_hits', 'BLAST_hits'])
        super()._check_input()

    def __set_informs(self) -> None:
        """
        Sets the informs.
        :return: None
        """
        self._informs['Min. percent identity'] = f"{self._parameters['min_percent_identity'].value}%"
        self._informs['Min. coverage'] = f"{self._parameters['min_coverage'].value}%"
        self._informs['_command'] = 'Gene detection mixed filtering'

    def __combine_hits(self) -> list[GeneDetectionMixedHit]:
        """
        Combines KMA + BLAST hits.
        :return: GeneDetectionMixedHit list
        """
        return [GeneDetectionMixedHit(self._tool_inputs['KMA_hits'][0].value.locus,
                                      self._tool_inputs['KMA_hits'][0].value.accession,
                                      self._tool_inputs['BLAST_hits'][0].value,
                                      self._tool_inputs['KMA_hits'][0].value)]

    def __create_summary_tables(self, combined_hits: list[GeneDetectionMixedHit]) -> list[pd.DataFrame]:
        """
        Creates a summary table based on genes discovered by both methods.
        :return: pandas dataframe
        """
        res_summary, res_complete = {}, {}
        score = 0
        for hit in combined_hits:
            the_locus = hit.locus
            the_accession = hit.accession
            kma_infos, blast_infos = [], []
            locus_length = 0
            if the_locus not in res_summary:
                res_summary[the_locus] = []
                res_complete[the_locus] = []
            if hit.kma_hit():
                kma_hit = hit.kma_hit()
                locus_length = kma_hit._length
                score += 1
                kma_infos = [kma_hit.percent_identity, kma_hit.subject_coverage, kma_hit._depth]
            if hit.blast_hit():
                blast_hit = hit.blast_hit()
                locus_length = blast_hit.blast_stats.subject_length
                score += 1
                blast_infos = [blast_hit.blast_stats.percent_identity, blast_hit.blast_stats.length_statistic]
            res_summary[the_locus] = [locus_length, 'str(the_accession)', score]
            res_complete[the_locus].extend(kma_infos)
            res_complete[the_locus].extend(blast_infos)
        return [pd.DataFrame.from_dict(res_summary, orient='index'),
                pd.DataFrame.from_dict(res_complete, orient='index')]
