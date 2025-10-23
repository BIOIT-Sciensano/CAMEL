import abc
import hashlib
from typing import Optional, Any, Union

from Bio.Seq import Seq

from camel.app.core.reports.htmlreportsection import HtmlReportSection


class TypingHitBase(metaclass=abc.ABCMeta):
    """
    This class represents is the base class for sequence typing hits. All hits should define the locus to which they
    belong and the allele that was detected. Abstract methods should be implemented by sub-classes so other classes
    can rely on these methods to provide the required functionality regardless of detection method.
    """

    SYMBOL_NO_HIT = '-'

    def __init__(self, locus: str, allele_id: str, new_allele_sequence: Optional[Seq] = None) -> None:
        """
        Initializes the typing hit.
        :param locus: Locus
        :param allele_id: Allele id of the hit
        :param new_allele_sequence: Sequence of the new allele
        """
        self._locus = locus
        self._allele_id = allele_id
        self._allele_page_url_template = None
        self._new_allele_sequence = new_allele_sequence

    @property
    def locus(self) -> str:
        """
        Returns the hit locus.
        :return: Locus
        """
        return self._locus

    @property
    def allele_id(self) -> str:
        """
        Returns the allele id.
        :return: Allele id
        """
        return self._allele_id

    @property
    def new_allele_hash(self) -> Union[str, None]:
        """
        Returns the hash of the novel allele (if available).
        :return: Hash (if available)
        """
        if not self.is_new_allele():
            return None
        return hashlib.md5(str(self._new_allele_sequence).encode('ascii')).hexdigest()[:6]

    def is_new_allele(self) -> bool:
        """
        Returns true if the allele is new.
        :return: True if new allele, False otherwise
        """
        return self._new_allele_sequence is not None

    @property
    def new_allele_sequence(self) -> Seq:
        """
        Returns the sequence of the novel allele sequence.
        """
        if not self.is_new_allele():
            raise ValueError('This hit is not a novel allele')
        return self._new_allele_sequence

    def set_allele_page_url_template(self, url_template: str) -> None:
        """
        Sets the allele page URL template.
        :param url_template: : URL template
        :return: None
        """
        self._allele_page_url_template = url_template

    @property
    def allele_page_url(self) -> Optional[str]:
        """
        Returns the allele info page URL if there is one.
        :return: Allele page URL
        """
        if self._allele_page_url_template is None:
            return None
        elif self._allele_id in ('?', '-'):
            return None
        return self._allele_page_url_template.format(allele_id=self.allele_id)

    @staticmethod
    @abc.abstractmethod
    def table_column_names() -> list[str]:
        """
        Returns the table column names.
        :return: Table column names
        """
        pass

    @abc.abstractmethod
    def to_table_row(self, hash_allele_ids: bool = False) -> list[str]:
        """
        Returns the hit as a row in a table.
        :param hash_allele_ids: If True, hashes for new allele ids are included
        :return: Table row
        """
        pass

    def to_dict(self, include_hashing: bool = False) -> dict[str, Any]:
        """
        Returns the hit as a dictionary.
        :param include_hashing: Whether to include the hashing info if a new allele is found
        :return: Hit as a dictionary
        """
        return {k: v for k, v in zip(self.table_column_names(), self.to_table_row())}

    @staticmethod
    @abc.abstractmethod
    def html_column_names() -> list[str]:
        """
        Returns the HTML column names.
        :return: HTML column names
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def to_html_row(self, report_section: HtmlReportSection, sub_dir: str = None) -> list[Any]:
        """
        Returns the hit as a row in a table.
        :param report_section: Section is passed to save the alignments
        :param sub_dir: Specific subdirectory of the base directory to store report files
        :return: Table row
        """
        pass

    @abc.abstractmethod
    def is_perfect_hit(self) -> bool:
        """
        Function to check if this is a perfect hit.
        :return: True if perfect hit, False otherwise
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def is_full_length(self) -> bool:
        """
        Function to check if this is a full length hit.
        :return: True if full length, False otherwise
        """
        raise NotImplementedError()

    @property
    @abc.abstractmethod
    def color(self) -> str:
        """
        Returns the color for this hit.
        :return: Color
        """
        raise NotImplementedError()
