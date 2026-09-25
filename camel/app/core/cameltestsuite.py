import logging
import os
import shutil
import tempfile
import time
import unittest
from pathlib import Path

from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.utils import reportutils

from camel.app.config import config
from camel.app.core.tool import Tool
from camel.app.loggers import logger


class CamelTestSuite(unittest.TestCase):
    """
    Base class for Camel test suites.
    """

    running_dir = None

    @staticmethod
    def get_test_file_dir(*args) -> Path:
        """
        Retrieves the directory with test files.
        :param args: Sub-directory or directories
        :return: Test files directory
        """
        dir_test = Path(config.dir_testdata, *args)
        if not dir_test.exists() or not dir_test.is_dir():
            raise FileNotFoundError(f"Cannot find test file directory: {dir_test}")
        return dir_test

    @staticmethod
    def get_reference_file_dir(*args) -> Path:
        """
        Retrieves the directory with reference files.
        :param args: Sub-directory or directories
        :return: Reference files directory
        """
        dir_ref = Path(config.dir_db, 'refgenomes', *args)
        if not dir_ref.exists() or not dir_ref.is_dir():
            raise FileNotFoundError(f"Cannot find reference file directory: {dir_ref}")
        return dir_ref

    def setUp(self) -> None:
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = Path(tempfile.mkdtemp(None, 'camel_', config.dir_temp))
        logger.debug(f"Directory for testing: {self.running_dir}")

    def tearDown(self) -> None:
        """
        Removes the working directory after executing the tests.
        :return: None
        """
        if os.environ.get('CAMEL_KEEP_TEST_DIRS') == '1':
            logger.debug(
                f"Keeping working directory (CAMEL_KEEP_TEST_DIRS): {self.running_dir}"
            )
            return
        if Path(self.running_dir).exists():
            logger.debug(f"Removing working directory: {self.running_dir}")
            CamelTestSuite._close_file_handlers(Path(self.running_dir))
            CamelTestSuite._remove_dir(Path(self.running_dir))

    @staticmethod
    def _close_file_handlers(dir_: Path) -> None:
        """
        Closes and detaches the logging file handlers that write into the given directory.
        Open file handles prevent removing the directory on NFS (which renames files before removing).
        :param dir_: Directory
        :return: None
        """
        dir_resolved = dir_.resolve()
        all_loggers = [logging.getLogger()] + [
            obj for obj in logging.Logger.manager.loggerDict.values() if isinstance(obj, logging.Logger)
        ]
        for logger_obj in all_loggers:
            for handler in list(logger_obj.handlers):
                if not isinstance(handler, logging.FileHandler):
                    continue
                if not Path(handler.baseFilename).resolve().is_relative_to(dir_resolved):
                    continue
                logger.warning(f"Closing log handler left open in test directory: {handler.baseFilename}")
                handler.close()
                logger_obj.removeHandler(handler)

    @staticmethod
    def _remove_dir(dir_: Path, attempts: int = 5, delay: float = 1.0) -> None:
        """
        Removes the given directory, retrying to give NFS time to release files that were just closed.
        :param dir_: Directory
        :param attempts: Number of attempts
        :param delay: Delay between attempts (in seconds)
        :return: None
        """
        for attempt in range(1, attempts + 1):
            try:
                shutil.rmtree(dir_)
                return
            except OSError as err:
                remaining = [str(p.relative_to(dir_)) for p in dir_.rglob('*') if p.is_file()] if dir_.exists() else []
                logger.warning(
                    f"Cannot remove directory (attempt {attempt}/{attempts}): {err}, remaining files: {remaining[:20]}"
                )
                if attempt == attempts:
                    raise
                time.sleep(delay)

    def verify_output_files(self, tool: Tool, key: str, nb_files: int = 1) -> None:
        """
        Verifies if the specified number of output files with the given key are created correctly.
        :param tool: Tool
        :param key: Output key
        :param nb_files: Number of output files that should be generated
        :return: None
        """
        self.assertIn(key, tool.tool_outputs, f"Key '{key}' missing from tool outputs")
        self.assertEqual(
            nb_files,
            len(tool.tool_outputs[key]),
            f"Unexpected number of tools outputs found: {len(tool.tool_outputs[key])}, expected {nb_files}",
        )
        for i in range(0, nb_files):
            output_file_path = tool.tool_outputs[key][i].path
            self.assertTrue(
                output_file_path.exists(),
                f"Output file '{key}' (index: {i}) does not exist",
            )
            self.assertGreater(
                output_file_path.stat().st_size,
                0,
                f"Output file '{key}' (index: {i}) is empty",
            )

    @staticmethod
    def export_report_section(report_section: HtmlReportSection, dir_out: Path) -> Path:
        """
        Exports the report section to the output directory.
        This method can be used to inspect the report generated by tools.
        :param dir_out: Output directory
        :param report_section: HTML report section
        :return: Path to the exported report
        """
        path_out = dir_out / 'report_section.html'
        path_out.parent.mkdir(exist_ok=True, parents=True)
        report = reportutils.init_report(
            path_out=path_out,
            key='Test report',
            title='Test report',
            dir_out=path_out.parent
        )
        report.add_html_object(report_section)
        report.save()
        logger.info(f'Report section exported to: {path_out}')
        return path_out
