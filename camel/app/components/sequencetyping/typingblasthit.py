from pathlib import Path
from typing import Optional, Any

from camel.app.components.blast.blasthitstatistics import BlastHitStatistics
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.components.sequencetyping.typinghitbase import TypingHitBase


class TypingBlastHit(TypingHitBase):
    """
    Sequence tying hit detected by blast.
    """

    SYMBOL_MULTI_HIT = '?'

    def __init__(self, locus: str, allele_id: str, type_, blast_stats: Optional[BlastHitStatistics]) -> None:
        """
        Initializes the hit.
        :param locus: Locus
        :param allele_id: Allele id
        :param type_: Locus type ('DNA', 'peptide')
        :param blast_stats: Blast hit statistics (optional)
        """
        novel_allele_seq = blast_stats.novel_allele_seq() if (
                blast_stats is not None and blast_stats.is_new_allele()) else None
        super().__init__(locus, allele_id, novel_allele_seq)
        self._type = type_
        self._blast_stats = blast_stats
        self._alignment_path = None

    @staticmethod
    def table_column_names() -> list[str]:
        """
        Returns the column names for the tabular output.
        :return: Table column names
        """
        return ['Locus', 'Allele', '% Identity', 'HSP/Locus length', 'Type']

    def to_table_row(self, hash_allele_ids: bool = False) -> list[str]:
        """
        Returns the hit as a row in a table.
        :param hash_allele_ids: If True, hashes for new allele ids are included
        :return: Table row
        """
        # Determine the allele id
        if not self.is_new_allele():
            allele_id = self.allele_id
        elif hash_allele_ids:
            allele_id = self.new_allele_hash
        else:
            allele_id = f'{self.allele_id}*'

        # Determine the % identity
        if self.is_new_allele() and hash_allele_ids is True:
            perc_identity = 100.0
        else:
            perc_identity = self._blast_stats.percent_identity if self.blast_stats else None

        # Return the output data
        return [
            self.locus,
            allele_id,
            f'{perc_identity:.2f}' if perc_identity is not None else '-',
            self.blast_stats.length_statistic if self.blast_stats else '-',
            self._type
        ]

    def to_dict(self, include_hashing: bool = False) -> dict[str, Any]:
        """
        Returns the hit as a dictionary.
        :param include_hashing: Whether to include the hashing info if a new allele is found
        :return: Hit dictionary
        """
        result = super().to_dict()
        if include_hashing and self.is_new_allele():
            result.update({
                'New allele': True,
                'Allele (hash)': self.new_allele_hash(full_length=True),
                'Allele sequence': str(self.new_allele_sequence)
            })
        return result

    @staticmethod
    def html_column_names() -> list[str]:
        """
        Returns the HTML column names.
        :return: HTML column names
        """
        return TypingBlastHit.table_column_names() + ['Alignment']

    def to_html_row(self, report_section: HtmlReportSection, sub_dir: Path = None) -> list[Any]:
        """
        Returns the hit as a row in a table.
        :param report_section: Section is passed to save the alignments
        :param sub_dir: Specific subdirectory of the base directory to store report files
        :return: HTML row
        """
        if self._alignment_path is None:
            alignment_cell = '-'
        else:
            relative_path = sub_dir / 'alignments' / self._alignment_path.name
            report_section.add_file(self._alignment_path, relative_path)
            alignment_cell = HtmlTableCell('view', link=str(relative_path))
        return [
            self.locus,
            HtmlTableCell(
                self.allele_id + ('*' if self.is_new_allele() else ''), self.color, link=self.allele_page_url),
            f'{self.blast_stats.percent_identity:.2f}' if self.blast_stats else '-',
            self.blast_stats.length_statistic if self.blast_stats else '-',
            self._type,
            alignment_cell
        ]

    @staticmethod
    def create_empty_hit(locus: str, type_: str) -> 'TypingBlastHit':
        """
        Returns an empty hit.
        :param locus: Locus
        :param type_: Locus type
        :return: None
        """
        return TypingBlastHit(locus, TypingHitBase.SYMBOL_NO_HIT, type_, None)

    @staticmethod
    def create_multi_hit(locus: str, type_: str, blast_stats: BlastHitStatistics) -> 'TypingBlastHit':
        """
        Returns a multi hit.
        :param locus: Locus
        :param type_: Locus type
        :param blast_stats: BLAST hit statistics
        :return: None
        """
        return TypingBlastHit(locus, TypingBlastHit.SYMBOL_MULTI_HIT, type_, blast_stats)

    @property
    def blast_stats(self) -> BlastHitStatistics:
        """
        Returns the BLAST statistics object.
        :return: BLAST hit statistics
        """
        return self._blast_stats

    @property
    def alignment_path(self) -> Path:
        """
        Returns the path to the alignment file.
        :return: Alignment file
        """
        return self._alignment_path

    @alignment_path.setter
    def alignment_path(self, alignment_path: Path) -> None:
        """
        Sets the alignment path.
        :param alignment_path: Alignment path
        :return: None
        """
        self._alignment_path = alignment_path

    @property
    def color(self) -> str:
        """
        Returns the color for this hit.
        Green: Perfect hit
        Light green: Full length hit with one or more mismatches
        Grey: Non-full length hit
        Red: No-hit
        :return: Color
        """
        if self.allele_id == TypingBlastHit.SYMBOL_MULTI_HIT:
            return 'yellow'
        elif self.allele_id == TypingBlastHit.SYMBOL_NO_HIT:
            return 'red'
        elif self.is_perfect_hit():
            return 'green'
        elif self.is_full_length():
            return 'lightgreen'
        return 'grey'

    def __repr__(self) -> str:
        """
        Returns the internal representation.
        :return: Representation
        """
        return f"TypingBlastHit('{self.locus}', allele_id='{self.allele_id}')"

    def is_perfect_hit(self) -> bool:
        """
        Function to check if this is a perfect hit.
        :return: True if perfect hit, False otherwise
        """
        if self._blast_stats is None:
            return False
        return self._blast_stats.is_perfect_hit()

    def is_new_allele(self) -> bool:
        """
        Checks if this hit is potentially a novel allele of the locus in the database.
        :return: True if the allele is new, False otherwise
        """
        if self._blast_stats is None:
            return False
        return self._blast_stats.is_new_allele()

    def is_full_length(self) -> bool:
        """
        Function to check if this is a full length hit.
        :return: True if full length, False otherwise
        """
        if self._blast_stats is None:
            return False
        return self._blast_stats.is_full_length()
