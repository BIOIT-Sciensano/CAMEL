from pathlib import Path
from typing import Any, Optional

from camelcore.app.reports.htmlelement import HtmlElement


class TsvExporter:
    """
    Exports data in tab separated values (TSV) format.
    """

    @staticmethod
    def export(output_data: list[list[Any]], header: Optional[list[str]], filename: str | Path,
               drop_columns: list[int] = None) -> None:
        """
        Exports output data in TSV format.
        :param output_data: Output data
        :param filename: Filename of the output file
        :param header: Header
        :param drop_columns: Columns not included in the exported file
        :return: None
        """
        if drop_columns is None:
            drop_columns = []
        with open(filename, 'w') as output_file:
            if header is not None:
                filtered_header = [h.strip() for i, h in enumerate(header) if i not in drop_columns]
                output_file.write('\t'.join(filtered_header))
                output_file.write('\n')
            for row in output_data:
                filtered_row = [TsvExporter.__get_element_text(e) for i, e in enumerate(row) if i not in drop_columns]
                output_file.write('\t'.join(filtered_row))
                output_file.write('\n')

    @staticmethod
    def __get_element_text(element: Any) -> str:
        """
        Returns the text from the given element.
        :param element: Element
        :return: Text
        """
        if isinstance(element, HtmlElement):
            return element.text
        return str(element)
