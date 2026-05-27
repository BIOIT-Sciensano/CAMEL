from dataclasses import dataclass

from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.tool import Tool
from camel.app.core.utils import toolutils
from camel.app.toolkits.export.tsvexporter import TsvExporter
from camel.app.toolkits.genedetection.genedetectionmixedhit import GeneDetectionMixedHit


@dataclass
class AggregatedHit:
    """
    This class represents a variant from a VCF file.
    """
    locus: str
    locus_length: int | str = '-'
    kma_color: str = 'red'
    blast_color: str = 'red'
    kma_identity: float | str = '-'
    kma_coverage: float | str = '-'
    kma_depth: float | str = '-'
    blast_identity: float | str = '-'
    blast_coverage: float | str = '-'
    blast_sequence: float | str = '-'
    blast_coordinates: float | str = '-'

    @staticmethod
    def header() -> list[str]:
        """
        Returns the header for the table.
        :return: header
        """
        return ['Locus', 'Locus length', 'KMA %identity', 'KMA %coverage', 'KMA depth', 'Blast %identity',
                'Blast HSP/Locus length', 'Blast Sequence (read or contig)', 'Blast Coordinates']

    def to_html_list(self) -> list:
        """
        Returns the gene detection information as a list.
        :return: list
        """
        return [HtmlTableCell(self.locus), HtmlTableCell(self.locus_length),
                HtmlTableCell(self.kma_identity, self.kma_color), HtmlTableCell(self.kma_coverage, self.kma_color),
                HtmlTableCell(self.kma_depth, self.kma_color), HtmlTableCell(self.blast_identity, self.blast_color),
                HtmlTableCell(self.blast_coverage, self.blast_color)]

    def to_list(self) -> list:
        """
        Returns the gene detection as a list.
        :return: list
        """
        return [self.locus, self.locus_length, self.kma_identity, self.kma_coverage, self.kma_depth,
                self.blast_identity, self.blast_coverage, self.blast_sequence, self.blast_coordinates]


class MixedHitAggregating(Tool):
    """
    Performs mixed hit aggregating from KMA and Blast gene detection.
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
        hits_combined = self.__combine_hits()
        aggregated_summary_hits = self.__create_combined_hits_summary(hits_combined)
        self._tool_outputs['VAL_hits'] = [ToolIOValue(hit) for hit in hits_combined]
        self._tool_outputs['VAL_agg_hits'] = [ToolIOValue(hit) for hit in aggregated_summary_hits]
        output_tsv_summary = self.folder / 'aggregated_hits_summary.tsv'
        TsvExporter.export([hit.to_list() for hit in aggregated_summary_hits],
                           aggregated_summary_hits[0].header(), str(output_tsv_summary))
        self._tool_outputs['TSV_agg_hits'] = [ToolIOFile(output_tsv_summary)]
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
        self._informs['_kma'] = self._input_informs['kma']
        self._informs['_blast'] = self._input_informs['blast']
        self._informs['Min. percent identity'] = f"{self._parameters['min_percent_identity'].value}%"
        self._informs['Min. coverage'] = f"{self._parameters['min_coverage'].value}%"
        self._informs['_command'] = 'Gene detection mixed filtering'

    def __combine_hits(self) -> list[GeneDetectionMixedHit]:
        """
        Combines KMA + BLAST hits into a GeneDetectionMixedHit list.
        :return: GeneDetectionMixedHit list
        """
        if not any([x for x in [self._tool_inputs['KMA_hits'], self._tool_inputs['BLAST_hits']]]):
            return []
        hits_per_locus = {}
        for hit in self._tool_inputs['KMA_hits']:
            locus = hit.value.locus
            if locus not in hits_per_locus:
                hits_per_locus[locus] = {'kma': None, 'blast': None}
            hits_per_locus[locus]['kma'] = hit.value
        for hit in self._tool_inputs['BLAST_hits']:
            locus = hit.value.locus
            if locus not in hits_per_locus:
                hits_per_locus[locus] = {'kma': None, 'blast': None}
            hits_per_locus[locus]['blast'] = hit.value
        output_list = [
            GeneDetectionMixedHit(locus=locus, accession=None, blast_hit=hits_per_locus[locus]['blast'],
                                  kma_hit=hits_per_locus[locus]['kma']) for locus in hits_per_locus
        ]
        return output_list

    def __create_combined_hits_summary(self, combined_hits: list[GeneDetectionMixedHit]) -> list[AggregatedHit]:
        """
        Creates a summary table based on genes discovered by both methods.
        :param combined_hits: List of GeneDetectionMixedHit objects
        :return: None
        """
        res_summary = []
        for hit in combined_hits:
            the_locus = self._retrieve_locus_name(hit.locus)
            agg_hit = AggregatedHit(locus=the_locus)
            kma_hit = hit.kma_hit()
            blast_hit = hit.blast_hit()
            if kma_hit:
                agg_hit.kma_identity = kma_hit.percent_identity
                agg_hit.kma_coverage = kma_hit.subject_coverage
                agg_hit.kma_depth = kma_hit.depth
                agg_hit.kma_color = kma_hit.color
                agg_hit.locus_length = kma_hit.length
            if blast_hit:
                agg_hit.blast_identity = blast_hit.blast_stats.percent_identity
                agg_hit.blast_coverage = blast_hit.blast_stats.length_statistic
                agg_hit.blast_sequence = blast_hit.blast_stats.query_id
                agg_hit.blast_coordinates = f'{blast_hit.blast_stats.query_start}..{blast_hit.blast_stats.query_end}'
                agg_hit.blast_color = blast_hit.color
                if agg_hit.locus_length == '-':
                    agg_hit.locus_length = blast_hit.blast_stats.subject_length
            res_summary.append(agg_hit)
        return res_summary

    def _retrieve_locus_name(self, locus: str) -> str:
        """
        Retrieves the locus name.
        :return: locus name
        """
        locus_name = self._input_informs['db']['mapping'].get_metadata(locus, 'allele')
        return locus_name
