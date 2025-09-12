import unittest

import pandas as pd

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltableformatter import HtmlTableFormatter
from camel.app.components.testing.cameltestsuite import CamelTestSuite


class TestHtmlTableFormatter(CamelTestSuite):
    """
    Tests the HTML table formatter.
    """

    def test_table_no_colors(self) -> None:
        """
        Tests regular formatting of a table.
        :return: None
        """
        # Sample data
        data_in = pd.DataFrame({
            'a': [1, 2500, 3, 400, 56000],
            'b': ['a', 'b', 'c', 'd', 'e'],
            'c': [10.2, 5.0, 12.42, 1355.4, 3.13451]
        })
        format_dict = [
            {'key': 'a', 'title': 'Column A', 'fmt': HtmlTableFormatter.INT_FMT},
            {'key': 'b', 'title': 'Column B'},
            {'key': 'c', 'title': 'Column C', 'fmt': HtmlTableFormatter.FLOAT_FMT}
        ]

        # Create the report section
        section = HtmlReportSection('Sample data')
        section.add_table(
            HtmlTableFormatter.format_table_data(data_in, format_dict),
            [c['title'] for c in format_dict],
            [('class', 'data')]
        )

        # Save the report in the current directory
        CamelTestSuite.export_report_section(section, self.running_dir)

    def test_table_with_colors(self) -> None:
        """
        Tests regular formatting of a table with colored rows.
        :return: None
        """
        # Sample data
        data_in = pd.DataFrame({
            'a': [1, 2500, 3, 400, 56000],
            'b': ['a', 'b', 'c', 'd', 'e'],
            'c': [10.2, 5.0, 12.42, 1355.4, 3.13451],
            'color': ['green', 'green', 'grey', 'red', 'green']
        })
        format_dict = [
            {'key': 'a', 'title': 'Column A', 'fmt': HtmlTableFormatter.INT_FMT},
            {'key': 'b', 'title': 'Column B'},
            {'key': 'c', 'title': 'Column C', 'fmt': HtmlTableFormatter.FLOAT_FMT}
        ]

        # Create the report section
        section = HtmlReportSection('Sample data')
        section.add_table(
            HtmlTableFormatter.format_table_data(data_in, format_dict, use_colors=True),
            [c['title'] for c in format_dict],
            [('class', 'data')]
        )

        # Save the report in the current directory
        CamelTestSuite.export_report_section(section, self.running_dir)

    def test_table_missing_column(self) -> None:
        """
        Tests that the function raises an error when a column is missing.
        :return: None
        """
        # Sample data
        data_in = pd.DataFrame({
            'a': [1, 2500, 3, 400, 56000],
            'b': ['a', 'b', 'c', 'd', 'e'],
            'c': [10.2, 5.0, 12.42, 1355.4, 3.13451]
        })
        format_dict = [
            {'key': 'a', 'title': 'Column A', 'fmt': HtmlTableFormatter.INT_FMT},
            {'key': 'b', 'title': 'Column B'},
            {'key': 'c', 'title': 'Column C', 'fmt': HtmlTableFormatter.FLOAT_FMT},
            {'key': 'non_existing', 'title': 'Non existing'}
        ]

        # Create the report section
        section = HtmlReportSection('Sample data')
        with self.assertRaises(ValueError):
            section.add_table(
                HtmlTableFormatter.format_table_data(data_in, format_dict),
                [c['title'] for c in format_dict],
                [('class', 'data')]
            )

    def test_table_missing_color(self) -> None:
        """
        Tests that the function raises an error when colored rows are requested but not provided.
        :return: None
        """
        # Sample data
        data_in = pd.DataFrame({
            'a': [1, 2500, 3, 400, 56000],
            'b': ['a', 'b', 'c', 'd', 'e'],
            'c': [10.2, 5.0, 12.42, 1355.4, 3.13451]
        })
        format_dict = [
            {'key': 'a', 'title': 'Column A', 'fmt': HtmlTableFormatter.INT_FMT},
            {'key': 'b', 'title': 'Column B'},
            {'key': 'c', 'title': 'Column C', 'fmt': HtmlTableFormatter.FLOAT_FMT},
        ]

        # Create the report section
        section = HtmlReportSection('Sample data')
        with self.assertRaises(ValueError):
            section.add_table(
                HtmlTableFormatter.format_table_data(data_in, format_dict, use_colors=True),
                [c['title'] for c in format_dict],
                [('class', 'data')]
            )


if __name__ == '__main__':
    unittest.main()
