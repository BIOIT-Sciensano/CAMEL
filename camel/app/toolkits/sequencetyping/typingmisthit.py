import json
from pathlib import Path
from typing import Any

from Bio.Seq import Seq
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltablecell import HtmlTableCell

from camel.app.toolkits.sequencetyping.typinghitbase import TypingHitBase


class TypingMiSTHit(TypingHitBase):
    """
    Sequence tying hit detected by blast.
    """

    def __init__(self, locus: str, allele_id: str, tags: list[str], novel_seq: str | None, allele_results: list) -> None:
        """
        Initializes the hit.
        :param locus: Locus
        :param allele_id: Allele id
        :param tags: Tags
        :param novel_seq: Novel allele sequence
        :return: None
        """
        seq = Seq(novel_seq) if novel_seq is not None else None
        super().__init__(locus, allele_id, new_allele_sequence=seq)
        self._tags = tags
        self._allele_results = allele_results

    @staticmethod
    def parse_mist_json(path_json: Path) -> list['TypingMiSTHit']:
        """
        Parses the hits from a MiST JSON output file.
        :param path_json: Path to the JSON file
        :return: List of hits
        """
        with path_json.open() as handle:
            data = json.load(handle)
        hits = []
        for locus, row in data['alleles'].items():
            hits.append(TypingMiSTHit(
                locus=locus,
                allele_id=row['allele_str'],
                tags=row['tags'],
                novel_seq=row['allele_results'][0]['sequence'] if len(row['allele_results']) > 0 else None,
                allele_results=row['allele_results']
            ))
        return hits

    @staticmethod
    def table_column_names() -> list[str]:
        """
        Returns the column names for the tabular output.
        :return: Table column names
        """
        return ['Locus', 'Allele', 'Tag(s)', 'Sequence', 'Position', 'Strand']

    def to_table_row(self, hash_allele_ids: bool = False) -> list[str]:
        """
        Returns the hit as a row in a table.
        :param hash_allele_ids: If True, hashes for new allele ids are included
        :return: Table row
        """
        return [
            self.locus,
            self.allele_id,
            ', '.join(self._tags) if len(self._tags) > 0 else '-',
            ' '.join(r['alignment']['seq_id'] for r in self._allele_results) if len(self._allele_results) > 0 else '-',
            ' '.join(f"{r['alignment']['start']}-{r['alignment']['end']}" for r in self._allele_results) if len(self._allele_results) > 0 else '-',
            ' '.join(r['alignment']['strand'] for r in self._allele_results)  if len(self._allele_results) > 0 else '-',
        ]

    @staticmethod
    def html_column_names() -> list[str]:
        """
        Returns the HTML column names.
        :return: HTML column names
        """
        return TypingMiSTHit.table_column_names()

    def to_html_row(self, report_section: HtmlReportSection, sub_dir: Path | None = None) -> list[Any]:
        """
        Returns the hit as a row in a table.
        :param report_section: Section is passed to save the alignments
        :param sub_dir: Specific subdirectory of the base directory to store report files
        :return: HTML row
        """
        return [
            self.locus,
            HtmlTableCell(self.allele_id, self.color, link=self.allele_page_url),
            ', '.join(self._tags) if len(self._tags) > 0 else '-',
            ' '.join(r['alignment']['seq_id'] for r in self._allele_results) if len(self._allele_results) > 0 else '-',
            ' '.join(f"{r['alignment']['start']}-{r['alignment']['end']}" for r in self._allele_results) if len(self._allele_results) > 0 else '-',
            ' '.join(r['alignment']['strand'] for r in self._allele_results) if len(self._allele_results) > 0 else '-',
        ]

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
        if self.is_perfect_hit():
            return 'green'
        return 'red'

    def __repr__(self) -> str:
        """
        Returns the internal representation.
        :return: Representation
        """
        return f"TypingMiSTHit('{self.locus}', allele_id='{self.allele_id}')"

    def is_perfect_hit(self) -> bool:
        """
        Function to check if this is a perfect hit.
        :return: True if perfect hit, False otherwise
        """
        return 'EXACT' in self._tags

    def is_full_length(self) -> bool:
        """
        Function to check if this is a full length hit.
        For MiST, full-length coverage is only confirmed for EXACT matches.
        :return: True if full length, False otherwise
        """
        return self.is_perfect_hit()
