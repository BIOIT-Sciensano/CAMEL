import json
from pathlib import Path
from typing import Any

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.components.sequencetyping.typinghitbase import TypingHitBase


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
        super().__init__(locus, allele_id, novel_seq)
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
                novel_seq=None, # TODO,
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
            ' '.join(r['alignment']['seq_id'] for r in self._allele_results),
            ' '.join(f"{r['alignment']['start']}-{r['alignment']['end']}" for r in self._allele_results),
            ' '.join(r['alignment']['strand'] for r in self._allele_results),
        ]

    def to_dict(self, include_hashing: bool = False) -> dict[str, Any]:
        """
        Returns the hit as a dictionary.
        :param include_hashing: Whether to include the hashing info if a new allele is found
        :return: Hit dictionary
        """
        result = super().to_dict()
        # if include_hashing and self.is_new_allele():
        #     result.update({
        #         'New allele': True,
        #         'Allele (hash)': self.new_allele_hash(full_length=True),
        #         'Allele sequence': str(self.new_allele_sequence)
        #     })
        return result

    @staticmethod
    def html_column_names() -> list[str]:
        """
        Returns the HTML column names.
        :return: HTML column names
        """
        return TypingMiSTHit.table_column_names()

    def to_html_row(self, report_section: HtmlReportSection, sub_dir: Path = None) -> list[Any]:
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
            ' '.join(r['alignment']['seq_id'] for r in self._allele_results),
            ' '.join(f"{r['alignment']['start']}-{r['alignment']['end']}" for r in self._allele_results),
            ' '.join(r['alignment']['strand'] for r in self._allele_results),
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
        return self._allele_id != TypingHitBase.SYMBOL_NO_HIT

    def is_new_allele(self) -> bool:
        """
        Checks if this hit is potentially a novel allele of the locus in the database.
        :return: True if the allele is new, False otherwise
        """
        return self._new_allele_sequence is not None

    def is_full_length(self) -> bool:
        """
        Function to check if this is a full length hit.
        :return: True if full length, False otherwise
        """
        if self.is_perfect_hit():
            return True
        return False
