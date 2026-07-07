from dataclasses import dataclass

from camelcore.app.reports.htmlelement import HtmlElement
from camelcore.app.reports.htmltablecell import HtmlTableCell


@dataclass(frozen=True)
class MixedHitHolder:
    """
    This class represents an aggregated hit with both KMA and Blast gene detection information.
    """
    cluster: str
    locus: str
    accession: str | HtmlTableCell
    locus_length: int | None = None
    kma_identity: float | None = None
    kma_coverage: float | None = None
    kma_depth: float | None = None
    blast_identity: float | None = None
    blast_coverage: str | None = None
    blast_sequence: str | None = None
    blast_coordinates: str | None = None
    kma_color: str = 'red'
    blast_color: str = 'red'

    @staticmethod
    def header_tsv() -> list[str]:
        """
        Returns the header for the TSV table.
        :return: header
        """
        return [
            'Cluster', 'Locus', 'Accession', 'Locus length', 'KMA %identity', 'KMA %coverage', 'KMA depth',
            'Blast %identity', 'Blast HSP/Locus length', 'Blast Sequence (read or contig)', 'Blast Coordinates'
        ]

    @staticmethod
    def header_html() -> list[str]:
        """
        Returns the header for the HTML table.
        :return: header
        """
        return [
            'Cluster', 'Locus', 'Accession', 'Locus length', '%identity', '%coverage', 'depth', '%identity',
            'HSP/Locus length'
        ]

    def kma_cells(self) -> list[HtmlTableCell | HtmlElement]:
        """
        Determines the type of cell to include in the report table.
        If there is a KMA hit, returns cells of %identity, coverage and depth.
        Alternatively, returns an element spanning three columns saying the hit was not detected by KMA.
        :return: List of HTML cells or a single HTML element
        """
        if self.kma_identity:
            return [
                HtmlTableCell(f'{self.kma_identity:.2f}', self.kma_color),
                HtmlTableCell(f'{self.kma_coverage:.2f}', self.kma_color),
                HtmlTableCell(f'{self.kma_depth:.2f}', self.kma_color)
            ]
        else:
            return [HtmlElement('td', 'Not detected by KMA.', [('colspan', 3)])]

    def blast_cells(self) -> list[HtmlTableCell | HtmlElement]:
        """
        Determines the type of Blast cell to include in the report table.
        If there is a blast hit, returns cells of %identity and coverage.
        Alternatively, returns an element spanning two columns saying the hit was not detected by Blast.
        :return: List of HTML cells or a single HTML element
        """
        if self.blast_identity:
            return [
                HtmlTableCell(f'{self.blast_identity:.2f}', self.blast_color),
                HtmlTableCell(self.blast_coverage, self.blast_color)
            ]
        else:
            return [HtmlElement('td', 'Not detected by BLAST.', [('colspan', 2)])]

    def to_html_row(self) -> list[HtmlTableCell | str]:
        """
        Returns the genes detected as an HTML list.
        :return: list
        """
        return [
            self.cluster,
            self.locus,
            self.accession,
            HtmlTableCell(self.locus_length),
        ] + self.kma_cells() + self.blast_cells()

    def to_table_row(self) -> list:
        """
        Returns the genes detected as a list.
        :return: list
        """
        return [
            self.cluster, self.locus, self.locus_length, self.accession, self.kma_identity, self.kma_coverage,
            self.kma_depth, self.blast_identity, self.blast_coverage, self.blast_sequence, self.blast_coordinates
        ]
