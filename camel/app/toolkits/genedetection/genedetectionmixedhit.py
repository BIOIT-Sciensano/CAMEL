from pathlib import Path
from typing import Union

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
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
        self.alignment_path = None
        self._cluster = None

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

    def is_perfect_hit(self) -> list:
        """
        Function to check if this is a perfect hit.
        :return: True if perfect hit, False otherwise
        """
        conditions = {'kma': self._kma_hit.is_perfect_hit() if self._kma_hit else False,
                      'blast': self._blast_hit.is_perfect_hit() if self._blast_hit else False}
        return conditions

    def is_full_length(self) -> list:
        """
        Function to check if this is a full length hit.
        :return: True if full length, False otherwise
        """
        conditions = {'kma': self._kma_hit.is_full_length() if self._kma_hit else False,
                      'blast': self._blast_hit.is_full_length() if self._blast_hit else False}
        return conditions

    @property
    def color(self) -> dict:
        """
        Returns the color for this hit.
        Green: Perfect hit
        Light green: Full length hit with one or more mismatches
        Grey: Non-full length hit
        :return: Color
        """
        colors = {'kma': self._kma_hit.color if self._kma_hit else 'red',
                  'blast': self._blast_hit.color if self._blast_hit else 'red'}
        return None
        # perfect_hits = self.is_perfect_hit()
        # full_length = self.is_full_length()
        # if perfect_hits[0] and perfect_hits[1]:
        #     return 'green'
        # elif full_length[0] and full_length[1]:
        #     return 'lightgreen'
        # else:
        #     return 'grey'

    @property
    def table_column_names(self) -> list[str]:
        """
        Returns the names of the columns of the tabular output.
        :return: List of column names
        """
        columns = ['DB_cluster', 'Locus', 'Accession', 'Detected by KMA', 'Detected by Blast']
        for metadata in self._metadata:
            columns.insert(-1, metadata['name'])
        return columns

    def to_table_row(self) -> list[str]:
        """
        Returns the hit as a table row.
        :return: List of table cell values
        """
        data = [
            self.cluster(),
            self.locus,
            self._accession if self._accession is not None else '-',
            ' ',
            ' '
        ]
        for metadata in self._metadata:
            data.insert(-2, metadata['value'])
        return data

    @property
    def html_column_names(self) -> list[str]:
        """
        Returns the names of the columns of the HTML output.
        :return: List of column names
        """
        return self.table_column_names[1:]

    def to_html_row(self, report_section: HtmlReportSection, sub_directory: Path, colored: bool = True) \
            -> list[Union[str, HtmlTableCell]]:
        """
        Returns the hit as a HTML table row.
        :param report_section: Section is passed to save additional data
        :param sub_directory: Subdirectory to save the additional data
        :param colored: If True, the row is colored
        :return: List of table cell values
        """
        print(self.to_table_row()[1:-3])
        print(self._get_accession_cell())
        test = [HtmlTableCell(v) for v in self.to_table_row()][1:-3] + [self._get_accession_cell()] + [
            HtmlTableCell(self.to_table_row()[-2], self._kma_hit.color if self._kma_hit else 'red')] + [
                   HtmlTableCell(self.to_table_row()[-1], self._blast_hit.color if self._blast_hit else 'red')]
        print(test)
        # print(self.color['kma'])
        # print(self.color['blast'])
        # exit()
        # return [HtmlTableCell(v, self.color) for v in self.to_table_row()][1:-1] + [self._get_accession_cell()]
        return ([HtmlTableCell(v) for v in self.to_table_row()][1:-3] + [self._get_accession_cell()] +
                [HtmlTableCell(self.to_table_row()[-2], self._kma_hit.color if self._kma_hit else 'red')] +
                [HtmlTableCell(self.to_table_row()[-1], self._blast_hit.color if self._blast_hit else 'red')])

    def kma_hit(self) -> GeneDetectionKMAHit:
        """
        Returns the KMA hit.
        :return: GeneDetectionKMAHit
        """
        return self._kma_hit

    def blast_hit(self) -> GeneDetectionBlastHit:
        """
        Returns the BLAST hit.
        :return:
        """
        return self._blast_hit

    def accession(self) -> str:
        """
        Returns the accession.
        :return: accession
        """
        return self._accession
