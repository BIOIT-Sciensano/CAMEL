import requests
from bs4 import BeautifulSoup


class PubMLSTParsingError(RuntimeError):
    """
    Error that is raised when parsing PubMLST fails.
    """
    pass


class PubMLSTParser:
    """
    Parsers for PubMLST web pages.
    """

    @staticmethod
    def parse_linked_data(url: str) -> tuple[str, str]:
        """
        Retrieves the linked data from a PubMLST page.
        :param url: URL of the web page to parse
        :return: Linked data
        """
        try:
            response = requests.get(url)
        except requests.exceptions.ConnectionError as err:
            raise PubMLSTParsingError(err)
        html = BeautifulSoup(PubMLSTParser.__cleanup_html(response.text), 'html.parser')
        try:
            html_value = html.find(text='Linked data').findNext('dd')
            html_label = html.find(text='Linked data').findNext('dt')
            return PubMLSTParser.__cleanup_label(html_label.text), html_value.text.replace(' PubMLST isolates', '')
        except AttributeError:
            raise PubMLSTParsingError('Cannot find linked data')

    @staticmethod
    def __cleanup_html(html_code):
        """
        Cleans up HTML code.
        :param html_code: HTML code
        :return: Cleaned up HTML code
        """
        return html_code.replace('<dd><', '<dd>&lt;')

    @staticmethod
    def __cleanup_label(label):
        """
        Capitalizes the label and replaces the underscore by a space.
        :param label: Input label
        :return: Cleaned up label
        """
        return label.capitalize().replace('_', ' ')
