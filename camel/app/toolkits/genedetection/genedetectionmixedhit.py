from dataclasses import replace
from pathlib import Path
from typing import Union

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.toolkits.genedetection.mixedhitholder import MixedHitHolder
from camel.app.toolkits.genedetection.genedetectionblasthit import GeneDetectionBlastHit
from camel.app.toolkits.genedetection.genedetectionhitbase import GeneDetectionHitBase
from camel.app.toolkits.genedetection.genedetectionkmahit import GeneDetectionKMAHit


class GeneDetectionMixedHit(GeneDetectionHitBase):
    """
    This class represents a gene detection hit detected with BLAST and/or KMA.
    """

    def __init__(self, locus: str, accession: str | None, blast_hit: GeneDetectionBlastHit | None,
                 kma_hit: GeneDetectionKMAHit | None) -> None:
        """
        Initializes this hit.
        :param locus: Locus
        :param accession: Accession
        :param blast_hit: Blast hit
        :param kma_hit: KMA hit
        :return: None
        """
        super().__init__(locus, accession)
        self._blast_hit = blast_hit
        self._kma_hit = kma_hit
        self._cluster = None

    @property
    def cluster(self) -> str | None:
        """
        Define cluster based on input.
        :return: Cluster as string
        """
        if self._kma_hit:
            return self._kma_hit.cluster
        elif self._blast_hit:
            return self._blast_hit.blast_stats.subject_id.split('__')[1]
        else:
            return None

    @property
    def color(self) -> str:
        """
        For a mixed hit, I depend on KMA and/or Blast colors. I set this to white so that the "accession" cell
        does not get colored by the "gene_detection_get_column_names" rule in the snakefile.
        :return: White
        """
        return 'white'

    def is_perfect_hit(self) -> dict[str, bool]:
        """
        Function to check if this is a perfect hit.
        :return: True if perfect hit, False otherwise
        """
        conditions = {'kma': self._kma_hit.is_perfect_hit() if self._kma_hit else False,
                      'blast': self._blast_hit.is_perfect_hit() if self._blast_hit else False}
        return conditions

    def is_full_length(self) -> dict[str, bool]:
        """
        Function to check if this is a full length hit.
        :return: True if full length, False otherwise
        """
        conditions = {'kma': self._kma_hit.is_full_length() if self._kma_hit else False,
                      'blast': self._blast_hit.is_full_length() if self._blast_hit else False}
        return conditions

    @property
    def table_column_names(self) -> list[str]:
        """
        Returns the names of the columns of the tabular output.
        :return: List of column names
        """
        columns = ['DB_cluster', 'Locus', 'Accession']
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
            self._accession if self._accession is not None else '-',
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
        return self.table_column_names[1:] + ['KMA', 'Blast']

    def to_html_row(self, report_section: HtmlReportSection, sub_directory: Path, colored: bool = True) \
            -> list[Union[str, HtmlTableCell]]:
        """
        Returns the hit as an HTML table row.
        :param report_section: Section is passed to save additional data
        :param sub_directory: Subdirectory to save the additional data
        :param colored: If True, the row is colored
        :return: List of table cell values
        """
        return ([HtmlTableCell(v) for v in self.to_table_row()][1:-1] + [self._get_accession_cell()] +
                [HtmlTableCell('v' if self._kma_hit else 'x', self._kma_hit.color if self._kma_hit else 'red')] +
                [HtmlTableCell('v' if self._blast_hit else 'x', self._blast_hit.color if self._blast_hit else 'red')])

    @property
    def holder(self) -> MixedHitHolder:
        """
        Returns this hit as a MixedHitHolder, combining the KMA and Blast statistics into a single object.
        :return: MixedHitHolder
        """
        mixed_hit = MixedHitHolder(cluster=self.cluster, locus=self.locus)
        if self._kma_hit:
            mixed_hit = replace(
                mixed_hit,
                kma_identity=self._kma_hit.percent_identity,
                kma_coverage=self._kma_hit.subject_coverage,
                kma_depth=self._kma_hit.depth,
                kma_color=self._kma_hit.color,
                locus_length=self._kma_hit.length)
        if self._blast_hit:
            blast_stats = self._blast_hit.blast_stats
            mixed_hit = replace(
                mixed_hit,
                blast_identity=blast_stats.percent_identity,
                blast_coverage=blast_stats.length_statistic,
                blast_sequence=blast_stats.query_id,
                blast_coordinates=f'{blast_stats.query_start}..{blast_stats.query_end}',
                blast_color=self._blast_hit.color,
                locus_length=blast_stats.subject_length if mixed_hit.locus_length is None else mixed_hit.locus_length)
        return mixed_hit

    def kma_hit(self) -> GeneDetectionKMAHit:
        """
        Returns the KMA hit class.
        :return: GeneDetectionKMAHit
        """
        return self._kma_hit

    def blast_hit(self) -> GeneDetectionBlastHit:
        """
        Returns the BLAST hit class.
        :return: GeneDetectionBlastHit
        """
        return self._blast_hit
