from camelcore.app.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool
from camel.app.core import toolutils
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
        all_hits_combined: dict[str, list[GeneDetectionMixedHit]] = self.__combine_hits()
        hits_combined = all_hits_combined['complete']
        hits_combined_summary = all_hits_combined['summary_hits']
        for hit in hits_combined:
            hit.accession = self._retrieve_accession_number(hit.locus)
            hit.locus = self._retrieve_locus_name(hit.locus)

        # Setting the outputs
        self._tool_outputs['VAL_hits'] = [ToolIOValue(hit) for hit in hits_combined_summary]
        self._tool_outputs['VAL_hits_complete'] = [ToolIOValue(hit) for hit in hits_combined]
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
        self._informs['extra_info'] = (
            'In the case of the mixed gene detection strategy, both <b>read-based (KMA)</b> and '
            '<b>assembly-based (Blast)</b> gene detection methods are used. Since Blast can report hits on multiple '
            'closely related alleles, and KMA only single best hits, '
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
    def __extract_locus_name_for_table(cluster_dic: dict) -> dict[str, str]:
        """
        Given a cluster dictionary, returns the locus with both KMA and Blast Hit.
        Alternatively, if either of those is empty, returns the first locus in appearance.
        :param cluster_dic: Cluster-Locus dictionary
        :return: List of locus/loci name(s)
        """
        res_locus = {'kma': [loc for loc in cluster_dic.keys() if cluster_dic[loc]['kma']],
                     'blast': [loc for loc in cluster_dic.keys() if cluster_dic[loc]['blast']]}

        # First: checking whether both methods detected the same locus in the cluster
        if res_locus['kma'] and res_locus['blast']:
            try:
                locus_to_return, extra_info = [loc for loc in res_locus['kma'] if loc in res_locus['blast']][0], 'both'
            except IndexError:
                # Returns one allele only but marking the extra info as diff - which will show in the HTML table
                # that both tools found different alleles
                locus_to_return, extra_info = res_locus['kma'][0], 'diff'
            return {'locus_to_return': locus_to_return, 'extra_info': extra_info}

        # Then: Checking which one detected a locus in the cluster.
        # Returns the first locus occurrence in the list.
        return {'locus_to_return': res_locus['kma'][0] if res_locus['kma'] else res_locus['blast'][0], 'extra_info': 'single'}

    def __build_hits_per_cluster_dictionary(self) -> dict[str, dict[str, dict[str, None]]]:
        """
        Builds the dictionary storing all loci hits by either tool, stored by cluster ID.
        That way, the method stores the unique hits of KMA vs potential multiple hits from Blast in the same cluster.
        :return: Dictionary
        """
        def __init_cluster_entry(
                dic_orig: dict,
                cluster_in: str,
                locus_in: str
        ) -> dict[str, dict[str, dict[str, None]]]:
            """
            Initialize the cluster entry in the dictionary if necessary.
            :param dic_orig: original dictionary
            :param cluster_in: cluster string
            :param locus_in: locus string
            :return: Initialized cluster dictionary or original dictionary
            """
            if cluster_in not in dic_orig:
                dic_orig[cluster_in] = {}
            if locus_in not in dic_orig[cluster_in].keys():
                dic_orig[cluster_in][locus_in] = {'kma': None, 'blast': None}
            return dic_orig

        hits_per_cluster = {}
        for hit in self._tool_inputs['HITS_kma']:
            cluster, locus = hit.value.cluster, hit.value.locus
            hits_per_cluster = __init_cluster_entry(
                dic_orig=hits_per_cluster,
                cluster_in=cluster,
                locus_in=locus
            )
            hits_per_cluster[cluster][locus]['kma'] = hit.value
        for hit in self._tool_inputs['HITS_blast']:
            blast_hit = hit.value
            updated_blast_hit = self.__set_blast_identity_same_as_kma(blast_hit)
            if updated_blast_hit is not None:
                cluster, locus = hit.value.cluster, hit.value.locus
                hits_per_cluster = __init_cluster_entry(
                    dic_orig=hits_per_cluster,
                    cluster_in=cluster,
                    locus_in=locus
                )
                hits_per_cluster[cluster][locus]['blast'] = hit.value
        return hits_per_cluster

    def __combine_hits(self) -> dict[str, list[GeneDetectionMixedHit]]:
        """
        Combines KMA + BLAST hits into a GeneDetectionMixedHit object list. The list is split in two for summary hits
        and the complete hit list.
        :return: GeneDetectionMixedHit dictionary
        """
        dictionary_hits_per_cluster = self.__build_hits_per_cluster_dictionary()
        output_summary_hits = []
        for cluster in dictionary_hits_per_cluster:
            loci_summary = self.__extract_locus_name_for_table(dictionary_hits_per_cluster[cluster])
            locus = loci_summary['locus_to_return']
            output_summary_hits.append(
                GeneDetectionMixedHit(
                    locus=locus,
                    accession=None,
                    blast_hit=dictionary_hits_per_cluster[cluster][locus]['blast'],
                    kma_hit=dictionary_hits_per_cluster[cluster][locus]['kma'],
                    extra_info=loci_summary['extra_info']
                )
            )
        all_loci_dictionary = {
            locus: dictionary_hits_per_cluster[cluster][locus]
            for cluster in dictionary_hits_per_cluster.keys()
            for locus in dictionary_hits_per_cluster[cluster].keys()
        }
        complete_hit_list = [
            GeneDetectionMixedHit(
                locus=locus,
                accession=None,
                blast_hit=all_loci_dictionary[locus]['blast'],
                kma_hit=all_loci_dictionary[locus]['kma']
            ) for locus in all_loci_dictionary
        ]
        return {'summary_hits': output_summary_hits, 'complete': complete_hit_list}

    def _retrieve_locus_name(self, locus: str) -> str:
        """
        Retrieves the locus name.
        :param locus: Locus
        :return: allele name
        """
        locus_name = self._input_informs['db']['mapping'].get_metadata(locus, 'allele')
        return locus_name

    def _retrieve_accession_number(self, locus: str) -> str:
        """
        Retrieves the accession number.
        :param locus: Locus
        :return: Accession cell
        """
        accession = self._input_informs['db']['mapping'].get_metadata(locus, 'accession', '-')
        return accession
