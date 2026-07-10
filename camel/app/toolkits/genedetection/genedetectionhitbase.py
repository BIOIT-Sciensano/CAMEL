import abc
from pathlib import Path
from typing import Any

from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltablecell import HtmlTableCell

from camel.app.toolkits.genedetection.genedetectionutils import GeneDetectionUtils
from camel.app.toolkits.genedetection.mapping import Mapping
from camel.app.toolkits.genedetection.mappingjson import MappingJSON


class GeneDetectionHitBase(metaclass=abc.ABCMeta):
    """
    This is the base class for hits detected by the gene detection workflows.
    """

    def __init__(self, locus: str, accession: str | None = None) -> None:
        """
        Initializes the hit.
        :param locus: Locus corresponding to the hit.
        :param accession: Accession number
        """
        self._locus = locus
        self._accession = accession
        self._metadata = []

    @property
    def locus(self) -> str:
        """
        Returns the locus.
        :return: Locus name
        """
        return self._locus

    @locus.setter
    def locus(self, locus: str) -> None:
        """
        Sets the locus for this hit.
        :param locus: Locus
        :return: None
        """
        self._locus = locus

    @property
    def accession(self) -> str | None:
        """
        Returns the accession.
        :return: Accession
        """
        return self._accession

    @accession.setter
    def accession(self, accession: str) -> None:
        """
        Setter for the accession.
        :param accession: Accession number
        :return: None
        """
        self._accession = accession

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
    def color(self) -> str:
        """
        Returns the color for this hit.
        Green: Perfect hit
        Light green: Full length hit with one or more mismatches
        Grey: Non-full length hit
        :return: Color
        """
        if self.is_perfect_hit():
            return 'green'
        elif self.is_full_length():
            return 'lightgreen'
        else:
            return 'grey'

    @property
    @abc.abstractmethod
    def table_column_names(self) -> list[str]:
        """
        Returns the names of the columns of the tabular output.
        :return: List of column names
        """
        raise NotImplementedError()

    @property
    def cluster(self) -> str:
        """
        Returns the cluster of the alignment.
        :return: Cluster
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def to_table_row(self) -> list[str]:
        """
        Returns the hit as a table row.
        :return: List of table cell values
        """
        raise NotImplementedError()

    def to_dict(self) -> dict[str, Any]:
        """
        Returns the hit as a dictionary.
        :return: Dictionary with hit data
        """
        return {k: v for k, v in zip(self.table_column_names, self.to_table_row())}

    @property
    @abc.abstractmethod
    def html_column_names(self) -> list[str]:
        """
        Returns the names of the columns of the HTML output.
        :return: List of column names
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def to_html_row(
        self, report_section: HtmlReportSection, sub_directory: Path, colored: bool = True
    ) -> list[str | HtmlTableCell]:
        """
        Returns the hit as a HTML table row.
        :param report_section: Section is passed to save additional data
        :param sub_directory: Subdirectory to save the additional data
        :param colored: If True, the row is colored
        :return: List of table cell values
        """
        raise NotImplementedError()

    def add_metadata(self, name: str, value: str) -> None:
        """
        Adds metadata to the hit.
        :param name: Metadata title
        :param value: Metadata value
        :return: None
        """
        self._metadata.append({'name': name, 'value': value})

    def resolve_names(self, mapping: Mapping | MappingJSON, extra_metadata: dict[str, str] | None = None) -> None:
        """
        Resolves this hit's locus name and accession from their standardized seq_id via the given mapping.
        :param mapping: Mapping of standardized seq_ids to original metadata
        :param extra_metadata: Optional extra per-DB metadata to add, as {'name': ..., 'key': ...}
        :return: None
        """
        seq_id = self.locus
        self.locus = mapping.get_metadata(seq_id, 'allele')
        self.accession = mapping.get_metadata(seq_id, 'accession', '-')
        if extra_metadata is not None:
            self.add_metadata(extra_metadata['name'], mapping.get_metadata(seq_id, extra_metadata['key']))

    def _get_accession_cell(self, colored: bool = True) -> HtmlTableCell:
        """
        Returns the table cell for the accession.
        :param colored: If True, colors the cell using this hit's color
        :return: Table cell with accession
        """
        color = self.color if colored else None
        if self._accession is None:
            return HtmlTableCell('-', color)
        elif GeneDetectionUtils.is_ncbi_accession(self._accession):
            link = f'https://www.ncbi.nlm.nih.gov/nuccore/{self._accession}'
            return HtmlTableCell(self._accession, color, link=link)
        else:
            return HtmlTableCell(self._accession, color)

    def get_metadata_value(self, key) -> Any:
        """
        Returns the value for the metadata with the given key.
        :param key: Key
        :return: Value
        """
        for dict_ in self._metadata:
            if dict_['name'] == key:
                return dict_['value']
        raise KeyError(key)
