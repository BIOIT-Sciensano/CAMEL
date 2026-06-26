from dataclasses import dataclass

from camel.app.core.reports.htmltablecell import HtmlTableCell


@dataclass(frozen=True)
class MixedHitHolder:
    """
    This class represents an aggregated hit with both KMA and Blast gene detection information.
    """
    cluster: str
    locus: str
    locus_length: int = None
    kma_identity: float = None
    kma_coverage: float = None
    kma_depth: float = None
    blast_identity: float = None
    blast_coverage: str = None
    blast_sequence: str = None
    blast_coordinates: str = None
    kma_color: str = 'red'
    blast_color: str = 'red'

    @staticmethod
    def header() -> list[str]:
        """
        Returns the header for the table.
        :return: header
        """
        return ['Cluster', 'Locus', 'Locus length', 'KMA %identity', 'KMA %coverage', 'KMA depth', 'Blast %identity',
                'Blast HSP/Locus length', 'Blast Sequence (read or contig)', 'Blast Coordinates']

    def to_html_row(self) -> list[HtmlTableCell | str]:
        """
        Returns the genes detected as an HTML list.
        :return: list
        """
        return [
            self.cluster,
            self.locus,
            HtmlTableCell(self.locus_length),
            HtmlTableCell(f'{self.kma_identity:.2f}' if self.kma_identity else '-', self.kma_color),
            HtmlTableCell(f'{self.kma_coverage:.2f}' if self.kma_coverage else '-', self.kma_color),
            HtmlTableCell(f'{self.kma_depth:.2f}' if self.kma_depth else '-', self.kma_color),
            HtmlTableCell(f'{self.blast_identity:.2f}' if self.blast_identity else '-', self.blast_color),
            HtmlTableCell(self.blast_coverage if self.blast_coverage else '-', self.blast_color)
        ]

    def to_table_row(self) -> list:
        """
        Returns the genes detected as a list.
        :return: list
        """
        return [
            self.cluster, self.locus, self.locus_length, self.kma_identity, self.kma_coverage, self.kma_depth,
            self.blast_identity, self.blast_coverage, self.blast_sequence, self.blast_coordinates
        ]
