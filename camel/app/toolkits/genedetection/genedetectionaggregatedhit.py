from dataclasses import dataclass

from camel.app.core.reports.htmltablecell import HtmlTableCell


@dataclass(frozen=True)
class AggregatedHit:
    """
    This class represents an aggregated hit with both KMA and Blast gene detection information.
    """
    cluster: str
    locus: str
    locus_length: int | str
    kma_identity: float | str
    kma_coverage: float | str
    kma_depth: float | str
    blast_identity: float | str
    blast_coverage: float | str
    blast_sequence: float | str
    blast_coordinates: float | str
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

    def to_html_list(self) -> list:
        """
        Returns the genes detected as an HTML list.
        :return: list
        """
        return [
            HtmlTableCell(self.cluster),
            HtmlTableCell(self.locus),
            HtmlTableCell(f'{self.locus_length:}'),
            HtmlTableCell(f'{self.kma_identity:.2f}' if isinstance(self.kma_identity, float) else self.kma_identity, self.kma_color),
            HtmlTableCell(f'{self.kma_coverage:.2f}' if isinstance(self.kma_coverage, float) else self.kma_coverage, self.kma_color),
            HtmlTableCell(f'{self.kma_depth:.2f}' if isinstance(self.kma_depth, float) else self.kma_depth, self.kma_color),
            HtmlTableCell(f'{self.blast_identity:.2f}' if isinstance(self.blast_identity, float) else self.blast_identity, self.blast_color),
            HtmlTableCell(f'{self.blast_coverage}', self.blast_color)
        ]

    def to_list(self) -> list:
        """
        Returns the genes detected as a list.
        :return: list
        """
        return [
            self.cluster, self.locus, self.locus_length, self.kma_identity, self.kma_coverage, self.kma_depth,
            self.blast_identity, self.blast_coverage, self.blast_sequence, self.blast_coordinates
        ]
