import dataclasses

from camelcore.app.io.tooliovalue import ToolIOValue

from camel.app.core import toolutils
from camel.app.core.tool import Tool
from camel.app.toolkits.genedetection.genedetectionblasthit import GeneDetectionBlastHit
from camel.app.toolkits.genedetection.genedetectionkmahit import GeneDetectionKMAHit
from camel.app.toolkits.genedetection.genedetectionmixedhit import GeneDetectionMixedHit


@dataclasses.dataclass(frozen=True)
class ClusterHits:
    """
    Represents detects hits for a given cluster.
    """

    cluster: str
    hits_kma: list[GeneDetectionKMAHit] = dataclasses.field(default_factory=list)
    hits_blast: list[GeneDetectionBlastHit] = dataclasses.field(default_factory=list)


class MixedHitAggregating(Tool):
    """
    Performs the aggregation of hits from KMA and Blast gene detection.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('Gene Detection: Mixed Hit Aggregating', '0.1')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Parse the input (creates lists of KMA and BLAST hits for each cluster)
        min_percent_identity = int(self.get_param_value('min_percent_identity'))
        min_coverage = int(self.get_param_value('min_coverage'))

        hits_by_cluster: dict[str, ClusterHits] = {}
        for hit in [io.value for io in self._tool_inputs['HITS_kma']]:
            if not hit.passes_thresholds(min_percent_identity, min_coverage):
                continue
            if hit.cluster not in hits_by_cluster:
                hits_by_cluster[hit.cluster] = ClusterHits(hit.cluster)
            hits_by_cluster[hit.cluster].hits_kma.append(hit)

        for hit in [io.value for io in self._tool_inputs['HITS_blast']]:
            if not hit.passes_thresholds(min_percent_identity, min_coverage):
                continue
            if hit.cluster not in hits_by_cluster:
                hits_by_cluster[hit.cluster] = ClusterHits(hit.cluster)
            hits_by_cluster[hit.cluster].hits_blast.append(hit)

        # Combines the KMA & BLAST hits into a single list of hits
        hits_summary: list[GeneDetectionMixedHit] = []
        for cluster, cluster_hits in hits_by_cluster.items():
            hits_summary.append(
                GeneDetectionMixedHit(
                    cluster=cluster, hits_blast=cluster_hits.hits_blast, hits_kma=cluster_hits.hits_kma
                )
            )
        self._tool_outputs['VAL_hits'] = [ToolIOValue(hit) for hit in hits_summary]
        self.__set_informs()

    def _check_input(self) -> None:
        """
        Checks whether the input is correct.
        :return: None
        """
        toolutils.check_input(self, keys_required=['HITS_kma', 'HITS_blast'])
        super()._check_input()

    def __set_informs(self) -> None:
        """
        Sets the informs.
        :return: None
        """
        self._informs['Min. percent identity'] = f"{int(self.get_param_value('min_percent_identity'))}%"
        self._informs['Min. coverage'] = f"{int(self.get_param_value('min_coverage'))}%"
        self._informs['extra_info'] = (
            'In the case of the mixed gene detection strategy, both <b>read-based (KMA)</b> and '
            '<b>assembly-based (BLAST)</b> gene detection methods are used.<ul>'
            '<li>BLAST can report hits on multiple closely related alleles, while KMA reports only a single '
            'best hit per cluster. The first table above shows one row per cluster, indicating whether it was '
            'detected by KMA, BLAST, or both (and whether they agreed on the same allele). The second table '
            'below shows the full detail for every allele detected by either method.</li>'
            '<li>BLAST computes percent identity over the aligned region only, while KMA normalizes over the '
            'full length of the reference sequence. To make the two directly comparable, BLAST\'s % identity and %'
            'covered values shown here are recalculated using the same principle as KMA (and capped at 100%).</li>'
            '</ul>'
        )
