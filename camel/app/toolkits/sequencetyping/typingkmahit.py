from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.toolkits.sequencetyping.typinghitbase import TypingHitBase


class TypingKMAHit(TypingHitBase):
    """
    Sequence tying hit detected by KMA.
    """

    def __init__(self, locus: str, allele_id: str, length: int, p_ident: float, p_cov: float, depth: float,
                 score: int) -> None:
        """
        Initializes the hit.
        :param locus: Locus
        :param allele_id: Allele id
        :param length: (Subject) length
        :param p_ident: (Subject) percent identity
        :param p_cov: (Subject) percent coverage
        :param depth: mean read depth across allele
        :param score: k-mer score
        """
        super().__init__(locus, allele_id)
        self._length = length
        self._p_ident = p_ident
        self._p_cov = p_cov
        self._depth = depth
        self._score = score

    @staticmethod
    def create_empty_hit(locus: str) -> 'TypingKMAHit':
        """
        Creates an empty KMA hit.
        :param locus: Locus
        :return: KMA hit
        """
        return TypingKMAHit(locus, '-', 0, 0, 0, 0, 0)

    @staticmethod
    def table_column_names() -> list[str]:
        """
        Returns the column names for the tabular output.
        :return: Table column names
        """
        return ['Locus', 'Allele', 'Length', '% Identity', '% Coverage', 'Depth']

    def to_table_row(self, hash_allele_ids: bool = False) -> list[str]:
        """
        Returns the hit as a row in a table.
        :param hash_allele_ids: If True, hashes for new allele ids are included
        :return: Table row
        """
        return [
            self.locus,
            self.allele_id,
            str(self._length) if self._length > 0 else '-',
            f'{float(self._p_ident):.2f}',
            f'{float(self._p_cov):.2f}',
            f'{float(self._depth):.2f}'
        ]

    @staticmethod
    def html_column_names() -> list[str]:
        """
        Returns the HTML column names.
        :return: HTML column names
        """
        return TypingKMAHit.table_column_names()

    def to_html_row(self, base_dir=None, sub_dir=None):
        """
        Returns the hit as a HTML row.
        :param base_dir: Base directory to store report
        :param sub_dir: Specific subdirectory of the base directory to store report files
        :return: Table row
        """
        return [
            self.locus,
            HtmlTableCell(self.allele_id, self.color, link=self.allele_page_url),
            str(self._length),
            f'{float(self._p_ident):.2f}',
            f'{float(self._p_cov):.2f}',
            f'{float(self._depth):.2f}',
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
        elif self.is_full_length():
            return 'lightgreen'
        else:
            return 'grey'

    def is_perfect_hit(self) -> bool:
        """
        Returns true if this is a perfect hit.
        :return: True if perfect
        """
        return self._p_ident == 100.0 and self.is_full_length()

    def is_full_length(self) -> bool:
        """
        Returns true if this is a full length hit.
        :return: True if full length
        """
        return self._p_cov == 100.0

    def __repr__(self) -> str:
        """
        Returns the internal representation.
        :return: Representation
        """
        return f"TypingKMAHit('{self.locus}', allele_id='{self.allele_id}')"

    @property
    def score(self) -> int:
        """
        Returns the k-mer score.
        :return: Score
        """
        return self._score

    @property
    def percent_identity(self) -> float:
        """
        Returns the percent identity of the hit.
        :return: Percent identity
        """
        return self._p_ident

    @property
    def percent_covered(self) -> float:
        """
        Returns the percentage of the subject sequence that is covered.
        :return:
        """
        return self._p_cov
