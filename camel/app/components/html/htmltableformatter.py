from typing import Union, TypedDict, Callable, Optional

import pandas as pd

from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.loggers import logger


class FormatEntry(TypedDict):
    """
    Helper class for type annotation of the formatting dictionary.
    """
    title: str
    key: str
    fmt: Optional[Callable]

class HtmlTableFormatter:
    """
    This class can be used to format Pandas dataframes as tables for the HTML output report.
    """
    FLOAT_FMT = lambda x: f'{x:.2f}'
    INT_FMT = lambda x: f'{int(x):,}'

    @staticmethod
    def _check_columns(data_in: pd.DataFrame, format_dict: list[FormatEntry], use_colors: bool) -> bool:
        """
        Checks if the specified columns are present in the input dataframe.
        :param data_in: Input dataframe
        :param format_dict: List of formatting dictionaries for each column
        :param use_colors: If True, the colors from the 'col' column are used for the rows
        :return: True if columns are present
        """
        for column in format_dict:
            if column['key'] in data_in.columns:
                continue
            raise ValueError(f"Expected column '{column['key']}' not found in dataframe.")
        if use_colors:
            if 'color' not in data_in.columns:
                raise ValueError("Expected column 'color' not found in dataframe.")
        return True

    @staticmethod
    def format_table_data(data_in: pd.DataFrame, format_dict: list[FormatEntry], use_colors: bool=False) -> \
            list[list[Union[HtmlTableCell, str]]]:
        """
        Formats the table data.
        :param data_in: Input dataframe
        :param format_dict: List of formatting dictionaries for each column
        :param use_colors: If True, the colors from the 'col' column are used for the rows
        :return: List of formatted table rows
        """
        logger.debug(f'Formatting table data ({len(data_in):,} rows)')
        HtmlTableFormatter._check_columns(data_in, format_dict, use_colors)
        table_data = []
        for row_in in data_in.to_dict('records'):
            current_row = []
            for column in format_dict:
                if column.get('fmt') is None:
                    current_row.append(row_in[column['key']])
                else:
                    try:
                        current_row.append(column['fmt'](row_in[column['key']]))
                    except ValueError as err:
                        logger.error(f'Failed to format {column["key"]}: {err}')
                        raise err
            if use_colors:
                row_color = row_in['color']
                current_row = [HtmlTableCell(value, color=row_color) for value in current_row]
            table_data.append(current_row)
        return table_data
