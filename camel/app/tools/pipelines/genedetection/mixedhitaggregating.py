from dataclasses import fields, replace

from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool
from camel.app.core.utils import toolutils
from camel.app.toolkits.export.tsvexporter import TsvExporter
from camel.app.toolkits.genedetection.genedetectionaggregatedhit import AggregatedHit
from camel.app.toolkits.genedetection.genedetectionblasthit import GeneDetectionBlastHit
from camel.app.toolkits.genedetection.genedetectionmixedhit import GeneDetectionMixedHit


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
        # Gathering and aggregating hits
        all_hits_combined = self.__combine_hits()
        hits_combined = all_hits_combined['complete']
        hits_combined_filt = all_hits_combined['unique_hits']
        aggregated_summary_hits = self.__create_combined_hits_summary(hits_combined)

        # Setting the outputs
        self._tool_outputs['VAL_hits'] = [ToolIOValue(hit) for hit in hits_combined_filt]
        self._tool_outputs['VAL_agg_hits'] = [ToolIOValue(hit) for hit in aggregated_summary_hits]
        output_tsv_summary = self.folder / 'aggregated_hits_summary.tsv'
        TsvExporter.export([hit.to_list() for hit in sorted(aggregated_summary_hits, key=lambda x: x.locus)],
                           aggregated_summary_hits[0].header(), str(output_tsv_summary))
        self._tool_outputs['TSV_agg_hits'] = [ToolIOFile(output_tsv_summary)]
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
        self._informs['Min. percent identity'] = f"{self._parameters['min_percent_identity'].value}%"
        self._informs['Min. coverage'] = f"{self._parameters['min_coverage'].value}%"
        self._informs['_extra_info'] = (
            'In the case of the mixed gene detection strategy, both <b>read-based (KMA)</b> and <b>assembly-based (Blast)</b> '
            'gene detection methods are used. Since Blast can report hits on multiple closely related alleles, and KMA only single best hits, '
            'the table contains matching alleles detected with both KMA and Blast. '
            'The details about all detected genes (i.e., all alleles) are presented in the second table.'
        )

    def __set_blast_identity_same_as_kma(self, hit: GeneDetectionBlastHit) -> GeneDetectionBlastHit | None:
        """
        Calculates and updates the blast identity to follow the same principle as KMA.
        If the new identity is below the identity threshold, the hit is now excluded further.
        :param hit: Blast hit
        :return: updated Blast hit
        """
        identity = hit.blast_stats.percent_identity
        subject_length = hit.blast_stats.subject_length
        subject_alignment_length = hit.blast_stats.alignment_length
        if subject_alignment_length != subject_length:
            new_identity = (subject_alignment_length / 100 * identity) / subject_length * 100
            if new_identity < float(self._parameters['min_percent_identity'].value):
                return None
            hit.blast_stats.percent_identity = new_identity
        return hit

    @staticmethod
    def __extract_locus_name_for_table(cluster_dic: dict) -> str:
        """
        Given a cluster list, returns the locus with both KMA and Blast Hit.
        Alternatively, if either of those is empty, returns the first locus in appearance.
        :return: Locus name
        """
        res_locus = {'kma': [loc for loc in cluster_dic.keys() if cluster_dic[loc]['kma']],
                     'blast': [loc for loc in cluster_dic.keys() if cluster_dic[loc]['blast']]}

        # First: checking whether both methods detected the same locus in the cluster
        if res_locus['kma'] and res_locus['blast']:
            return [loc for loc in res_locus['kma'] if loc in res_locus['blast']][0]

        # Then: Checking which one detected a locus in the cluster.
        # Returns the first locus occurrence in the list.
        return res_locus['kma'][0] if res_locus['kma'] else res_locus['blast'][0]

    def __combine_hits(self) -> dict[str, list[GeneDetectionMixedHit]]:
        """
        Combines KMA + BLAST hits into a GeneDetectionMixedHit object list.
        :return: GeneDetectionMixedHit list
        """
        # If no hits (this shouldn't happen), an exception is thrown
        if not any([x for x in [self._tool_inputs['HITS_kma'], self._tool_inputs['HITS_blast']]]):
            raise ValueError("Empty hit tables. This shouldn't happen.")

        def __init_cluster_entry(dic_orig: dict, cluster_in: str, locus_in: str) -> dict[str, dict[str, dict[str, None]]]:
            """
            Initialize the cluster entry in the dictionary if necessary.
            :param dic_orig: original dictionary
            :param cluster_in: cluster string
            :param locus_in: locus string
            :return: initialized cluster dictionary or original dictionary
            """
            if cluster_in not in dic_orig:
                dic_orig[cluster] = {}
            if locus_in not in dic_orig[cluster_in].keys():
                dic_orig[cluster_in][locus_in] = {'kma': None, 'blast': None}
            return dic_orig

        # Now looping over all hits to store hits, arranged by clusters.
        # That way, it catches unique hits of KMA vs multiple hits from Blast in the same cluster.
        hits_per_cluster = {}
        for hit in self._tool_inputs['HITS_kma']:
            cluster, locus = hit.value.cluster, hit.value.locus
            hits_per_cluster = __init_cluster_entry(dic_orig=hits_per_cluster, cluster_in=cluster, locus_in=locus)
            hits_per_cluster[cluster][locus]['kma'] = hit.value
        for hit in self._tool_inputs['HITS_blast']:
            blast_hit = hit.value
            updated_blast_hit = self.__set_blast_identity_same_as_kma(blast_hit)
            if updated_blast_hit is not None:
                cluster, locus = hit.value.cluster, hit.value.locus
                hits_per_cluster = __init_cluster_entry(dic_orig=hits_per_cluster, cluster_in=cluster, locus_in=locus)
                hits_per_cluster[cluster][locus]['blast'] = hit.value
        filter_unique_hits_dic = {}
        for cluster in hits_per_cluster:
            locus_to_add = self.__extract_locus_name_for_table(hits_per_cluster[cluster])
            filter_unique_hits_dic[locus_to_add] = hits_per_cluster[cluster][locus_to_add]
        output_unique_hits = [
            GeneDetectionMixedHit(
                locus=locus,
                accession=None,
                blast_hit=filter_unique_hits_dic[locus]['blast'],
                kma_hit=filter_unique_hits_dic[locus]['kma']
            ) for locus in filter_unique_hits_dic
        ]
        all_loci_dictionary = {
            locus: hits_per_cluster[cluster][locus]
            for cluster in hits_per_cluster.keys()
            for locus in hits_per_cluster[cluster].keys()
        }
        complete_hit_list = [
            GeneDetectionMixedHit(
                locus=locus,
                accession=None,
                blast_hit=all_loci_dictionary[locus]['blast'],
                kma_hit=all_loci_dictionary[locus]['kma']
            ) for locus in all_loci_dictionary
        ]
        return {'unique_hits': output_unique_hits, 'complete': complete_hit_list}

    def __build_aggregated_hit(self, hit: GeneDetectionMixedHit) -> AggregatedHit:
        """
        Takes a mixed hit object and transform into an AggregatedHit object.
        :param hit: GeneDetectionMixedHit object
        :return: AggregatedHit object
        """
        agg_hit = AggregatedHit(**{f.name: '-' for f in fields(AggregatedHit) if
                                   f.name not in ['kma_color', 'blast_color']})
        agg_hit = replace(
            agg_hit,
            locus=self._retrieve_locus_name(hit.locus),
            cluster=hit.cluster()
        )
        if kma_hit := hit.kma_hit():
            agg_hit = replace(
                agg_hit,
                kma_identity=kma_hit.percent_identity,
                kma_coverage=kma_hit.subject_coverage,
                kma_depth=kma_hit.depth,
                kma_color=kma_hit.color,
                locus_length=kma_hit.length
            )
        if blast_hit := hit.blast_hit():
            blast_stats = blast_hit.blast_stats
            agg_hit = replace(
                agg_hit,
                blast_identity=f'{blast_stats.percent_identity:.2f}',
                blast_coverage=blast_stats.length_statistic,
                blast_sequence=blast_stats.query_id,
                blast_coordinates=f'{blast_stats.query_start}..{blast_stats.query_end}',
                blast_color=blast_hit.color,
                locus_length=agg_hit.locus_length if agg_hit.locus_length != '-' else blast_stats.subject_length
            )
        return agg_hit

    def __create_combined_hits_summary(self, combined_hits: list[GeneDetectionMixedHit]) -> list[AggregatedHit]:
        """
        Creates a list containing all aggregated hits.
        :param combined_hits: List of GeneDetectionMixedHit objects
        :return: List of all the aggregated hits
        """
        return [self.__build_aggregated_hit(hit) for hit in combined_hits]

    def _retrieve_locus_name(self, locus: str) -> str:
        """
        Retrieves the locus name.
        :param locus: Locus
        :return: allele name
        """
        locus_name = self._input_informs['db']['mapping'].get_metadata(locus, 'allele')
        return locus_name
