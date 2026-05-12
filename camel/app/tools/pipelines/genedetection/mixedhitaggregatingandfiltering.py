from pathlib import Path

import pandas as pd

from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool
from camel.app.core.utils import toolutils
from camel.app.toolkits.genedetection.genedetectionmixedhit import GeneDetectionMixedHit


class MixedHitAggregatingAndFiltering(Tool):
    """
    Performs mixed hit aggregating and filtering (KMA + BLAST).
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
        output_tsv_name = self.folder / 'aggregated_hits.tsv'
        self.__create_summary_table(hits_combined, output_tsv_name)
        self._tool_outputs['VAL_hits'] = [ToolIOValue(hit) for hit in hits_combined]
        self._tool_outputs['TSV_aggregated'] = [ToolIOFile(output_tsv_name)]
        self.__set_informs()

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
        if self._tool_inputs['KMA_hits'] == [] and self._tool_inputs['BLAST_hits'] == []:
            return []
        dictionary_hits_per_locus = {}
        for h in self._tool_inputs['KMA_hits']:
            locus = h.value.locus
            if locus not in dictionary_hits_per_locus:
                dictionary_hits_per_locus[locus] = {'kma': [], 'blast': []}
            dictionary_hits_per_locus[locus]['kma'].append(h)
        for h in self._tool_inputs['BLAST_hits']:
            locus = h.value.locus
            if locus not in dictionary_hits_per_locus:
                dictionary_hits_per_locus[locus] = {'kma': [], 'blast': []}
            dictionary_hits_per_locus[locus]['blast'].append(h)
        output_list = []
        for locus in dictionary_hits_per_locus:
            kma_hit, blast_hit, accession = None, None, None
            if dictionary_hits_per_locus[locus]['kma']:
                accession = dictionary_hits_per_locus[locus]['kma'][0].value.accession
                kma_hit = dictionary_hits_per_locus[locus]['kma'][0].value
            if dictionary_hits_per_locus[locus]['blast']:
                accession = accession or dictionary_hits_per_locus[locus]['blast'][0].value.accession
                blast_hit = dictionary_hits_per_locus[locus]['blast'][0].value
            output_list.append(GeneDetectionMixedHit(locus, accession, blast_hit, kma_hit))
        return output_list

    def __create_summary_table(self, combined_hits: list[GeneDetectionMixedHit], output_path: Path) -> list | None:
        """
        Creates a summary table based on genes discovered by both methods.
        :param combined_hits: List of GeneDetectionMixedHit objects
        :param output_path: Path to output directory
        :return: None
        """
        res_summary = []
        header_list = ['Locus', 'Locus length', 'KMA %identity', 'KMA %coverage', 'KMA depth', 'Blast %identity',
                       'Blast HSP/Locus length', 'Blast Sequence (read or contig)', 'Blast Coordinates', 'Color']
        if not combined_hits:
            return res_summary
        for hit in combined_hits:
            info_to_store = {'locus': '-', 'locus_length': '-', 'kma_id': '-', 'kma_cov': '-', 'kma_depth': '-',
                             'blast_id': '-', 'blast_cov': '-', 'blast_sequence': '-', 'blast_coordinates': '-', 'color': '-'}
            the_locus, the_accession = self._retrieve_locus_name(hit.locus)
            kma_hit = hit.kma_hit()
            blast_hit = hit.blast_hit()
            info_to_store['locus'] = the_locus
            info_to_store['kma_id'] = '-' if not kma_hit else kma_hit.percent_identity
            info_to_store['kma_cov'] = '-' if not kma_hit else kma_hit.subject_coverage
            info_to_store['kma_depth'] = '-' if not kma_hit else kma_hit.depth
            info_to_store['locus_length'] = blast_hit.blast_stats.subject_length if not kma_hit else kma_hit.length
            info_to_store['blast_id'] = '-' if not blast_hit else blast_hit.blast_stats.percent_identity
            info_to_store['blast_cov'] = '-' if not blast_hit else blast_hit.blast_stats.length_statistic
            info_to_store['blast_sequence'] = '-' if not blast_hit else blast_hit.blast_stats.query_id,
            info_to_store['blast_coordinates'] = '-' if not blast_hit else f'{blast_hit.blast_stats.query_start}..{blast_hit.blast_stats.query_end}',
            info_to_store['color'] = hit.color
            res_summary.append(info_to_store)
        summary_df = pd.DataFrame(sorted(res_summary, key=lambda x: x['locus']))
        summary_df.columns = header_list
        summary_df.to_csv(output_path, sep='\t', index=False)

    def _retrieve_locus_name(self, locus: str) -> tuple[str, str]:
        """
        Retrieves the locus names.
        :return: locus names
        """
        updated_hit = self._input_informs['db']['mapping'].get_metadata(locus, 'allele')
        updated_accession = self._input_informs['db']['mapping'].get_metadata(locus, 'accession', '-')
        return updated_hit, updated_accession
