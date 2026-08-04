import unittest

from camelcore.app.reports.htmlexpandablediv import HtmlExpandableDiv
from camelcore.app.reports.htmlreport import HtmlReport

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.loggers import logger


class TestHtmlReporter(CamelTestSuite):
    """
    Tests the HtmlReport and related classes.
    """

    def test_html_report(self) -> None:
        """
        Tests the creation and export of a simple HtmlReport.
        :return: None
        """
        report_path = self.running_dir / 'report.html'
        report = HtmlReport(report_path, self.running_dir)
        report.initialize("Test report")
        report.add_header("Test report", 1)
        report.add_paragraph("This  is the report content")
        report.save()
        logger.info(f"Report saved in: {report_path}")
        self.assertGreater(report_path.stat().st_size, 0)

    def test_report_with_js(self) -> None:
        """
        Tests the creation of a report that includes Javascript.
        :return: None
        """
        report_path = self.running_dir / 'report.html'
        report = HtmlReport(report_path, self.running_dir)
        report.initialize("Test report")
        report.add_header("Test report", 1)
        div = HtmlExpandableDiv('large-table', 'Table')
        table_data = [['row', i] for i in range(0, 40)]
        div.add_table(table_data, ['Col 1', 'Col 2'], [('class', 'data')])
        report.add_html_object(div)
        report.save()
        logger.info(f"Report saved in: {report_path}")
        self.assertGreater(report_path.stat().st_size, 0)

    def test_pipeline_header(self):
        """
        Tests adding a pipeline header to a report.
        :return: None
        """
        report_path = self.running_dir / 'report.html'
        report = HtmlReport(report_path, self.running_dir)
        report.initialize("Test report")
        report.add_pipeline_header('My <i>pipeline</i>')
        report.save()
        logger.info(f"Report saved in: {report_path}")
        self.assertGreater(report_path.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
