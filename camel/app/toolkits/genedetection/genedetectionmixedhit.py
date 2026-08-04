import itertools
from pathlib import Path

from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltablecell import HtmlTableCell

from camel.app.toolkits.genedetection.genedetectionblasthit import GeneDetectionBlastHit
from camel.app.toolkits.genedetection.genedetectionhitbase import GeneDetectionHitBase
from camel.app.toolkits.genedetection.genedetectionkmahit import GeneDetectionKMAHit
from camel.app.toolkits.genedetection.mapping import Mapping
from camel.app.toolkits.genedetection.mappingjson import MappingJSON


class GeneDetectionMixedHit(GeneDetectionHitBase):
    """
    This class represents a gene detection hit detected with BLAST and/or KMA.
    """

    def __init__(
        self,
        cluster: str,
        hits_blast: list[GeneDetectionBlastHit] | None,
        hits_kma: list[GeneDetectionKMAHit] | None,
    ) -> None:
        """
        Initializes this hit.
        :param hits_blast: List of blast hit
        :param hits_kma: List of KMA hits
        :return: None
        """
        locus = GeneDetectionMixedHit.extract_locus_name(hits_blast, hits_kma)
        super().__init__(locus, accession=None)
        self._blast_hits = hits_blast if hits_blast is not None else []
        self._kma_hits = hits_kma if hits_kma is not None else []
        self._cluster = cluster

        # Groups by locus (based on seq_id, before the names are mapped to the original loci)
        self._hits_by_locus = GeneDetectionMixedHit.__group_hits_by_locus(self._blast_hits, self._kma_hits)
        self._kma_matches_blast = bool(self._kma_hits) and self._hits_by_locus[self._kma_hits[0].locus]['blast'] is not None

    @staticmethod
    def extract_locus_name(hits_blast: list[GeneDetectionBlastHit], hits_kma: list[GeneDetectionKMAHit]) -> str:
        """
        Extracts the locus name.
        :param hits_blast: Detected BLAST hits
        :param hits_kma: Detected KMA hits
        :return: Locus name
        """
        hits: list[GeneDetectionHitBase] = hits_blast + hits_kma
        locus_names = set(h.locus for h in hits)
        return ', '.join(sorted(locus_names))

    @staticmethod
    def __group_hits_by_locus(
        hits_blast: list[GeneDetectionBlastHit], hits_kma: list[GeneDetectionKMAHit]
    ) -> dict[str, dict[str, GeneDetectionBlastHit | GeneDetectionKMAHit | None]]:
        """
        Groups the BLAST and KMA hits by the (standardized) seq_id of their locus, so that a KMA hit and a BLAST
        hit referring to the same underlying sequence can later be rendered on a single row.
        :param hits_blast: Detected BLAST hits
        :param hits_kma: Detected KMA hits
        :return: Dictionary of seq_id -> {'kma': hit or None, 'blast': hit or None}
        """
        hits_by_locus: dict[str, dict[str, GeneDetectionBlastHit | GeneDetectionKMAHit | None]] = {}
        for hit in hits_kma:
            hits_by_locus.setdefault(hit.locus, {'kma': None, 'blast': None})['kma'] = hit
        for hit in hits_blast:
            hits_by_locus.setdefault(hit.locus, {'kma': None, 'blast': None})['blast'] = hit
        return hits_by_locus

    def resolve_names(self, mapping: Mapping | MappingJSON, extra_metadata: dict[str, str] | None = None) -> None:
        """
        Resolves the locus name, accession and extra metadata of each underlying KMA/Blast hit via the given
        mapping, then rebuilds this hit's combined locus name (and extra metadata) from the resolved values.
        :param mapping: Mapping of standardized seq_ids to original metadata
        :param extra_metadata: Optional extra per-DB metadata to add, as {'name': ..., 'key': ...}
        :return: None
        """
        for hit in itertools.chain(self._kma_hits, self._blast_hits):
            hit.resolve_names(mapping, extra_metadata)
        self._locus = GeneDetectionMixedHit.extract_locus_name(self._blast_hits, self._kma_hits)
        if extra_metadata is not None:
            values = sorted({
                h.get_metadata_value(extra_metadata['name'])
                for h in itertools.chain(self._kma_hits, self._blast_hits)
            })
            self.add_metadata(extra_metadata['name'], ', '.join(values))

    @property
    def cluster(self) -> str | None:
        """
        Returns the DB cluster.
        :return: Cluster as string
        """
        return self._cluster

    @property
    def color(self) -> str:
        """
        Returns the color for the hit.
        :return: Color
        """
        if self._blast_hits and self._kma_hits:
            return 'green' if self._kma_matches_blast else 'lightgreen'
        if self._blast_hits or self._kma_hits:
            return 'yellow'
        return 'red'

    def is_perfect_hit(self) -> bool:
        """
        Function to check if this is a perfect hit.
        This is a placeholder not used as colors are set depending on the respective detection by either BLAST or KMA.
        :return: Dictionary of boolean
        """
        return all(h.is_perfect_hit() for h in itertools.chain(self._kma_hits, self._blast_hits))

    def is_full_length(self) -> bool:
        """
        Function to check if this is a full length hit.
        This is a placeholder not used as colors are set depending on the respective detection by either BLAST or KMA.
        :return: True if full length, False otherwise
        """
        return all(h.is_full_length() for h in itertools.chain(self._kma_hits, self._blast_hits))

    @property
    def detected_by(self) -> str:
        """
        Returns the detected by string.
        :return: Detected by string
        """
        if self._blast_hits and self._kma_hits:
            return 'Both' if self._kma_matches_blast else 'Both - diff. alleles'
        if self._blast_hits:
            return 'BLAST only'
        if self._kma_hits:
            return 'KMA only'
        return 'None'

    @property
    def table_column_names(self) -> list[str]:
        """
        Returns the columns names for the tabular output.
        :return: List of column names
        """
        columns = ['DB_cluster', 'Locus', 'Detected by']
        for metadata in self._metadata:
            columns.insert(-1, metadata['name'])
        return columns

    def to_table_row(self) -> list[str]:
        """
        Returns the hit as a table row.
        :return: List of table cell values
        """
        data = [
            self.cluster,
            self.locus,
            self.detected_by,
        ]
        for metadata in self._metadata:
            data.insert(-1, metadata['value'])
        return data

    @property
    def html_column_names(self) -> list[str]:
        """
        Returns the names of the columns of the HTML output.
        :return: List of column names
        """
        return self.table_column_names[1:-1] + ['Detected by']

    def to_html_row(
        self, report_section: HtmlReportSection, sub_directory: Path, colored: bool = True
    ) -> list[str | HtmlTableCell]:
        """
        Returns the hit as an HTML table row.
        :param report_section: Section is passed to save additional data
        :param sub_directory: Subdirectory to save the additional data
        :param colored: If True, the row is colored
        :return: List of table cell values
        """
        row = [self.locus, HtmlTableCell(self.detected_by, self.color)]
        for metadata in self._metadata:
            row.insert(-1, metadata['value'])
        return row

    @property
    def kma_hits(self) -> list[GeneDetectionKMAHit]:
        """
        Returns the KMA hit class (if available).
        :return: GeneDetectionKMAHit
        """
        return self._kma_hits

    @property
    def blast_hits(self) -> list[GeneDetectionBlastHit]:
        """
        Returns the BLAST hit class (if available).
        :return: GeneDetectionBlastHit
        """
        return self._blast_hits

    @property
    def hits_by_locus(self) -> dict[str, dict[str, GeneDetectionBlastHit | GeneDetectionKMAHit | None]]:
        """
        Returns the hits grouped by locus, combining the KMA and Blast hits that refer to the same sequence.
        :return: Dictionary of seq_id -> {'kma': hit or None, 'blast': hit or None}
        """
        return self._hits_by_locus
