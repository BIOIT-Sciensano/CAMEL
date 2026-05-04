from typing import Any

import pandas as pd

from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.reports.htmltablecell import HtmlTableCell
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
        Combines KMA + BLAST hits into a GeneDetectionMixedHit list.
        :return: GeneDetectionMixedHit list
        """
        return [GeneDetectionMixedHit(self._tool_inputs['KMA_hits'][0].value.locus,
                                      self._tool_inputs['KMA_hits'][0].value.accession,
                                      None if self._tool_inputs['BLAST_hits'] == [] else self._tool_inputs['BLAST_hits'][0].value,
                                      None if self._tool_inputs['KMA_hits'] == [] else self._tool_inputs['KMA_hits'][0].value)]

    def __create_summary_tables(self, combined_hits: list[GeneDetectionMixedHit]) \
            -> list[list[Any] | Any]:
        """
        Creates a summary table based on genes discovered by both methods.
        :return: pandas dataframe
        """
        res_summary, res_complete = [], {}
        res_summary_to_return = {'header': ['Locus', 'Locus length', 'Accession'], 'values': []}

        for hit in combined_hits:
            the_locus = hit.locus
            if the_locus not in res_complete:
                res_complete[the_locus] = []
            if hit.kma_hit():
                kma_hit = hit.kma_hit()
                locus_length = kma_hit.length
                res_complete[the_locus].extend([kma_hit.percent_identity, kma_hit.subject_coverage, kma_hit.depth])
            else:
                res_complete[the_locus].extend(['-', '-', '-'])
            if hit.blast_hit():
                blast_hit = hit.blast_hit()
                locus_length = blast_hit.blast_stats.subject_length
                res_complete[the_locus].extend([
                    blast_hit.blast_stats.percent_identity, blast_hit.blast_stats.length_statistic])
            else:
                res_complete[the_locus].extend(['-', '-'])
            res_summary = [the_locus, locus_length, hit.accession]
            res_summary_to_return['values'].append([HtmlTableCell(v, hit.color) for v in res_summary])
        return [res_summary_to_return,
                pd.DataFrame.from_dict(res_complete, orient='index')]
