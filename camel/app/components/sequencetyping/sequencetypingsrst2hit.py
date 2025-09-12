from typing import Any, Optional

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.components.sequencetyping.sequencetypinghitbase import SequenceTypingHitBase


class SequenceTypingSRST2Hit(SequenceTypingHitBase):
    """
    Sequence tying hit detected by SRST2.
    """

    def __init__(self, locus: str, allele_id: str, mismatches: str, uncertainty: str, depth: Optional[float]) -> None:
        """
        Initializes the hit.
        :param locus: Locus
        :param allele_id: Allele id
        :param mismatches: Mismatches between reads and allele
        :param uncertainty: Uncertainty between reads and allele
        :param depth: mean read depth across allele
        """
        super().__init__(locus, allele_id)
        self._mismatches = mismatches
        self._uncertainty = uncertainty
        self._depth = depth

    @staticmethod
    def create_empty_hit(locus: str) -> 'SequenceTypingSRST2Hit':
        """
        Creates an empty SRST2 hit.
        :param locus: Locus
        :return: SRST2 hit
        """
        return SequenceTypingSRST2Hit(locus, '-', '-', '-', 0.0)

    @staticmethod
    def table_column_names() -> list[str]:
        """
        Returns the column names for the tabular output.
        :return: Table column names
        """
        return ['Locus', 'Allele', 'Mismatches', 'Uncertainty', 'Depth']

    def to_table_row(self, hash_allele_ids: bool = False) -> list[str]:
        """
        Returns the hit as a row in a table.
        :param hash_allele_ids: If True, hashes for new allele ids are included
        :return: Table row
        """
        return [
            self.locus,
            self.allele_id,
            self._mismatches,
            self._uncertainty,
            f'{float(self._depth):.2f}' if self._depth is not None else '-'
        ]

    @staticmethod
    def html_column_names() -> list[str]:
        """
        Returns the HTML column names.
        :return: HTML column names
        """
        return SequenceTypingSRST2Hit.table_column_names()

    def to_html_row(self, report_section: HtmlReportSection, sub_dir: str = None) -> list[Any]:
        """
        Returns the hit as a HTML row.
        :param report_section: Section is passed to save the alignments
        :param sub_dir: Specific subdirectory of the base directory to store report files
        :return: HTML row
        """
        return [
            self.locus,
            HtmlTableCell(self.allele_id, self.color, link=self.allele_page_url),
            self._mismatches,
            self._uncertainty,
            f'{float(self._depth):.2f}' if self._depth is not None else '-'
        ]

    @property
    def color(self) -> str:
        """
        Returns the color for this hit.
        :return: Color
        """
        if self.allele_id == '-':
            return 'red'
        elif self.is_perfect_hit():
            return 'green'
        elif self._uncertainty == "-":
            return 'lightgreen'
        else:
            return 'grey'

    def is_perfect_hit(self) -> bool:
        """
        Function to check if this is a perfect hit.
        :return: True if perfect hit, False otherwise
        """
        return self._mismatches == "0" and self._uncertainty == "-"

    def is_full_length(self) -> bool:
        """
        Function to check if this is a full length hit.
        :return: True if full length, False otherwise
        """
        return 'hole' in self._mismatches

    def __repr__(self) -> str:
        """
        Returns the internal representation.
        :return: Representation
        """
        return f"TypingSRST2Hit('{self.locus}', allele_id='{self.allele_id}', mismatch='{self._mismatches}')"
