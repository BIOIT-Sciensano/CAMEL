from dataclasses import dataclass

from camel.app.core.reports.htmltablecell import HtmlTableCell


@dataclass
class AggregatedHit:
    """
    This class represents an aggregated hit with both KMA and Blast gene detection information.
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
        Returns the genes detected as an HTML list.
        :return: list
        """
        return [
            HtmlTableCell(self.locus),
            HtmlTableCell(f'{self.locus_length}:'),
            HtmlTableCell(f'{self.kma_identity}:.2f', self.kma_color),
            HtmlTableCell(f'{self.kma_coverage}:.2f', self.kma_color),
            HtmlTableCell(f'{self.kma_depth}:2f', self.kma_color),
            HtmlTableCell(f'{self.blast_identity}:2f', self.blast_color),
            HtmlTableCell(f'{self.blast_coverage}:.2f', self.blast_color)
        ]

    def to_list(self) -> list:
        """
        Returns the genes detected as a list.
        :return: list
        """
        return [
            self.locus, self.locus_length, self.kma_identity, self.kma_coverage, self.kma_depth,
            self.blast_identity, self.blast_coverage, self.blast_sequence, self.blast_coordinates
        ]
