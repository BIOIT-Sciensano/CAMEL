import unittest

import shutil

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
from camel.app.tools.samtools.samtoolsflagstat import SamtoolsFlagstat
from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
from camel.app.tools.samtools.samtoolsindexcram import SamtoolsIndexCram
from camel.app.tools.samtools.samtoolsmerge import SamtoolsMerge
from camel.app.tools.samtools.samtoolssort import SamtoolsSort
from camel.app.tools.samtools.samtoolsview import SamtoolsView


class TestSamtools(CamelTestSuite):
    """
    Tests the Samtools tool suite.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('samtools')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'reference.fasta')
    FILE_BAM = ToolIOFile(test_file_dir / 'aln.bam')
    FILE_BAM1 = ToolIOFile(test_file_dir / 'aln1.bam')
    FILE_BAM2 = ToolIOFile(test_file_dir / 'aln2.bam')
    FILE_CRAM = ToolIOFile(test_file_dir / 'aln.cram')
    FILE_TXTdepth = ToolIOFile(test_file_dir / 'aln_depth.txt')

    def test_samtools_depth(self) -> None:
        """
        Test SamtoolsDepth
        :return: None
        """
        samtools_depth = SamtoolsDepth()
        samtools_depth.add_input_files({'BAM': [TestSamtools.FILE_BAM]})
        samtools_depth.run(self.running_dir)
        self.verify_output_files(samtools_depth, 'TSV')
        self.assertIn('median_depth', samtools_depth.informs)

    def test_samtools_fastaindex(self) -> None:
        """
        Tests SamtoolsFastaIndex.
        :return: None
        """
        samtools_fastaindex = SamtoolsFastaIndex()
        samtools_fastaindex.add_input_files({'FASTA': [TestSamtools.FILE_FASTA_REF]})
        samtools_fastaindex.update_parameters(symlink_input=True)
        samtools_fastaindex.run(self.running_dir)
        self.verify_output_files(samtools_fastaindex, 'FASTA')
        self.assertTrue(samtools_fastaindex.tool_outputs['FASTA'][0].path.is_symlink())

    def test_samtools_fastaindex_no_symlink(self) -> None:
        """
        Tests SamtoolsFastaIndex.
        :return: None
        """
        path_in = self.running_dir / 'copied_ref.fasta'
        shutil.copyfile(TestSamtools.FILE_FASTA_REF.path, path_in)
        samtools_fastaindex = SamtoolsFastaIndex()
        samtools_fastaindex.add_input_files({'FASTA': [ToolIOFile(path_in)]})
        samtools_fastaindex.update_parameters(symlink_input=False)
        samtools_fastaindex.run(self.running_dir)
        self.verify_output_files(samtools_fastaindex, 'FASTA')
        self.assertFalse(samtools_fastaindex.tool_outputs['FASTA'][0].path.is_symlink())

    def test_samtools_flagstat(self) -> None:
        """
        Tests SamtoolsFlagstat.
        :return: None
        """
        samtools_flagstat = SamtoolsFlagstat()
        samtools_flagstat.add_input_files({'BAM': [TestSamtools.FILE_BAM]})
        samtools_flagstat.run(self.running_dir)
        self.verify_output_files(samtools_flagstat, 'TXT')
        informs_expected = [
            'total', 'secondary', 'supplementary', 'duplicates', 'mapped', 'paired', 'read1', 'read2',
            'properly_paired', 'singletons', 'mapping_rate']
        for inform in informs_expected:
            self.assertIn(inform, samtools_flagstat.informs, f"{inform} not found in informs.")

    def test_samtools_index(self) -> None:
        """
        Tests SamtoolsIndex.
        :return: None
        """
        samtools_index = SamtoolsIndex()
        samtools_index.add_input_files({'BAM': [TestSamtools.FILE_BAM]})
        samtools_index.run(self.running_dir)
        self.verify_output_files(samtools_index, 'BAM')

    def test_samtools_index_cram(self) -> None:
        """
        Tests SamtoolsIndexCram.
        :return: None
        """
        samtools_index_cram = SamtoolsIndexCram()
        samtools_index_cram.add_input_files({
            'CRAM': [TestSamtools.FILE_CRAM],
            'FASTA_REF': [TestSamtools.FILE_FASTA_REF]
        })
        samtools_index_cram.run(self.running_dir)
        self.verify_output_files(samtools_index_cram, 'CRAI')

    def test_samtools_merge(self) -> None:
        """
        Tests SamtoolsMerge.
        :return: None
        """
        samtools_merge = SamtoolsMerge()
        samtools_merge.add_input_files({'BAM': [TestSamtools.FILE_BAM1, TestSamtools.FILE_BAM2]})
        samtools_merge.run(self.running_dir)
        self.verify_output_files(samtools_merge, 'BAM')

    def test_samtools_sort(self) -> None:
        """
        Tests SamtoolsSort.
        :return: None
        """
        samtools_sort = SamtoolsSort()
        samtools_sort.add_input_files({'BAM': [TestSamtools.FILE_BAM]})
        samtools_sort.run(self.running_dir)
        self.verify_output_files(samtools_sort, 'BAM')

    def test_samtools_view(self) -> None:
        """
        Tests SamtoolsView.
        :return: None
        """
        samtools_view = SamtoolsView()
        samtools_view.add_input_files({'BAM': [TestSamtools.FILE_BAM]})
        samtools_view.update_parameters(output_format='SAM')
        samtools_view.run(self.running_dir)
        self.verify_output_files(samtools_view, 'SAM')


if __name__ == '__main__':
    unittest.main()
