import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.strainge.straingstkmerize import StrainGSTKmerize
from camel.app.tools.strainge.straingstreporter import StrainGSTReporter
from camel.app.tools.strainge.straingstrun import StrainGSTRun


class TestStrainGST(CamelTestSuite):
    """
    Initializes this testing tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('strainge')
    FILE_FASTQ = ToolIOFile(test_file_dir / 'bsubtilis.fastq.gz')
    FILE_HDF5 = ToolIOFile(test_file_dir / 'bsubtilis.hdf5')
    DB_HDF5 = ToolIOFile(test_file_dir / 'small_bacillus_db.hdf5')

    def test_straingst_kmerize(self) -> None:
        """
        Testing StrainGST kmerize on ont sequencing data.
        :return: None
        """
        straingst_kmerize = StrainGSTKmerize()
        straingst_kmerize.add_input_files({'FASTQ': [TestStrainGST.FILE_FASTQ]})
        straingst_kmerize.run(self.running_dir)
        self.verify_output_files(straingst_kmerize, 'HDF5')

    def test_straingst_run(self) -> None:
        """
        Testing StrainGST run on ONT sequencing data and a small B. subtilis database.
        :return: None
        """
        straingst_run = StrainGSTRun()
        straingst_run.add_input_files({'HDF5': [TestStrainGST.FILE_HDF5], 'DB_HDF5': [TestStrainGST.DB_HDF5]})
        straingst_run.run(self.running_dir)
        self.verify_output_files(straingst_run, 'TSV_STATS')

    def test_straingst_reporter(self) -> None:
        """
        Tests StrainGST followed by the reporter.
        :return: None
        """
        # First run straingst
        straingst_run = StrainGSTRun()
        straingst_run.add_input_files({'HDF5': [TestStrainGST.FILE_HDF5], 'DB_HDF5': [TestStrainGST.DB_HDF5]})
        straingst_run.run(self.running_dir)
        self.verify_output_files(straingst_run, 'TSV_STATS')

        # Then run the reporter
        reporter = StrainGSTReporter()
        reporter.add_input_files({'TSV': straingst_run.tool_outputs['TSV_STRAINS']})
        reporter.add_input_informs({'straingst': straingst_run.informs})
        reporter.update_parameters(suffix='ont')
        reporter.run(self.running_dir)
        output_section = reporter.tool_outputs['VAL_HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
