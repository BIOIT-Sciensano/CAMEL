import unittest
from pathlib import Path

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.utils import fastqutils


class TestFastqUtils(CamelTestSuite):
    """
    Tests the Fastq utils module.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('components')

    def test_get_sample_name_miseq_fmt(self) -> None:
        """
        Tests the get sample name function for MiSEQ format.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(Path('/data/temp/Z4686_S31_L001_R1_001.fastq.gz')), 'Z4686')

    def test_get_sample_name_miseq_fmt_no_s(self) -> None:
        """
        Tests the get sample name function for MiSEQ format without the sample number.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(Path('/data/temp/Z4686_L001_R1_001.fastq.gz')), 'Z4686')

    def test_get_sample_name_miseq_fmt_no_s_lowercase(self) -> None:
        """
        Tests the get sample name function for MiSEQ format without the sample number.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(Path('/data/temp/Z4686_L001_R1_001.fastq.gz')), 'Z4686')

    def test_get_sample_name_miseq_fmt_no_s_no_l(self) -> None:
        """
        Tests the get sample name function for MiSEQ format without the sample number and lane number.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(Path('/data/temp/Z4686_R1_001.fastq.gz')), 'Z4686')

    def test_get_sample_name_simple_fmt(self) -> None:
        """
        Tests the get sample name function for the simple format.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(Path('/data/temp/Z4686_1.fastq.gz')), 'Z4686')

    def test_get_sample_name_simple_fmt_alt(self) -> None:
        """
        Tests the get sample name function for the simple format.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(Path('/data/temp/S18BD02705_R1.fastq.gz')), 'S18BD02705')

    def test_get_sample_name_ont_bacdis(self) -> None:
        """
        Tests the get sample name function for the bacdis ONT format.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(
            Path('/data/temp/S22BD04543_SQK-RBK114-96_SUPv5-0-0_2025_run0113_B25.fastq.gz'), fastqutils.PATTERN_FQ_ONT), 'S22BD04543')
        self.assertEqual(fastqutils.get_sample_name(
            Path('/data/temp/S22BD04543_SQK-RPB114-24_HACv5-0-0_2025_run0113_B25.fastq'), fastqutils.PATTERN_FQ_ONT), 'S22BD04543')
        # Regular filename (not matching the format)
        self.assertEqual(fastqutils.get_sample_name(
            Path('/data/temp/S22BD04543_ont.fastq'), fastqutils.PATTERN_FQ_ONT),'S22BD04543_ont')
        self.assertEqual(fastqutils.get_sample_name(
            Path('/data/temp/Myco-SRR8948399_ont-ds.fastq.gz'), fastqutils.PATTERN_FQ_ONT),'Myco-SRR8948399_ont-ds')

    def test_get_sample_name_invalid_fmt(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        with self.assertRaises(ValueError):
            fastqutils.get_sample_name(Path('/data/temp/invalid_name.fastq.gz'))

    def test_get_sample_name_with_dots(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(Path('/data/temp/my-sample.1.fastq.gz')), 'my-sample')

    def test_get_sample_name_fq_ext(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(Path('/data/temp/my-sample_1.fq')), 'my-sample')

    def test_get_sample_name_fq_ext_gzipped(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(Path('/data/temp/my-sample_1.fq.gz')), 'my-sample')

    def test_get_sample_name_with_p(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(Path('/data/temp/my-sample_1P.fastq')), 'my-sample')

    def test_get_sample_name_with_underscore(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(Path('/data/temp/my_sample_1.fastq')), 'my_sample')

    def test_get_sample_name_se(self) -> None:
        """
        Tests the get sample name function for a single end sample name.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(
            Path('/data/temp/my_sample.fastq'), fastqutils.PATTERN_FQ_SE), 'my_sample')

    def test_get_sample_name_se_gzipped(self) -> None:
        """
        Tests the get sample name function for a single end sample name.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(
            Path('/data/temp/my_sample.fastq.gz'), fastqutils.PATTERN_FQ_SE), 'my_sample')

    def test_get_sample_name_se_gzipped_with_dashes(self) -> None:
        """
        Tests the get sample name function for a single end sample name.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(
            Path('/temp/my_sample22-ds.fastq.gz'), fastqutils.PATTERN_FQ_SE), 'my_sample22-ds')

    def test_get_sample_name_parentheses(self) -> None:
        """
        Tests the get sample name function that contains parentheses.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(Path('/data/temp/UI-586(SRR7648453)_1.fastq.gz')), 'UI-586SRR7648453')

    def test_get_sample_name_miseq(self) -> None:
        """
        Tests the get sample name function that contains parentheses.
        :return: None
        """
        self.assertEqual(fastqutils.get_sample_name(
            Path('/data/temp/S20BD03018_S16_L001_R1_001.fastq.gz')), 'S20BD03018')

    def test_count_reads_interleaved(self) -> None:
        """
        Tests function that returns the number of reads.
        :return: None
        """
        self.assertEqual(fastqutils.count_reads(TestFastqUtils.test_file_dir / 'interleaved.fq'), 4)

    def test_count_reads_interleaved_gzipped(self) -> None:
        """
        Tests function that returns the number of reads.
        :return: None
        """
        self.assertEqual(fastqutils.count_reads(TestFastqUtils.test_file_dir / 'interleaved.fq.gz'), 4)

    def test_get_all_read_names(self) -> None:
        """
        Test the function that retrieves the read names from the given fastq file.
        :return: None
        """
        read_names = fastqutils.get_all_read_names(TestFastqUtils.test_file_dir / 'fq-fwd.fq')
        self.assertEqual(read_names, {
            'M04115:7:000000000-AMA6W:1:1101:18190:1854', 'M04115:7:000000000-AMA6W:1:1101:23504:1916',
            'M04115:7:000000000-AMA6W:1:1101:23856:1955', 'M04115:7:000000000-AMA6W:1:1101:9914:1907'})

    def test_get_all_read_names_gzipped_input(self) -> None:
        """
        Test the function that retrieves the read names from the given gzipped fastq file.
        :return: None
        """
        read_names = fastqutils.get_all_read_names(TestFastqUtils.test_file_dir / 'fq-fwd.fq.gz')
        self.assertEqual(read_names, {
            'M04115:7:000000000-AMA6W:1:1101:18190:1854', 'M04115:7:000000000-AMA6W:1:1101:23504:1916',
            'M04115:7:000000000-AMA6W:1:1101:23856:1955', 'M04115:7:000000000-AMA6W:1:1101:9914:1907'})

    def test_count_bases(self) -> None:
        """
        Test the function to count the number of bases in a fastq file
        :return: None
        """
        count = fastqutils.count_bases(TestFastqUtils.test_file_dir / 'fq-pe.fq')
        self.assertEqual(count, 1502)

    def test_count_bases_gzip_input(self) -> None:
        """
        Test the function to count the number of bases in a gzipped fastq file
        :return: None
        """
        count = fastqutils.count_bases(TestFastqUtils.test_file_dir / 'fq-pe.fq.gz')
        self.assertEqual(count, 1502)

    def test_is_fastq_with_fastq(self) -> None:
        """
        Tests the function that checks whether the input file is a FASTQ file using a FASTQ file as input.
        :return: None
        """
        input_file = TestFastqUtils.test_file_dir / 'fq-file.fastq'
        is_fastq = fastqutils.is_fastq(input_file)
        self.assertTrue(is_fastq)

    def test_is_fastq_with_gzip_fastq(self) -> None:
        """
        Tests the function that checks whether the input file is a FASTQ file using a gzipped FASTQ file as input.
        :return: None
        """
        input_file = TestFastqUtils.test_file_dir / 'fq-file.fastq.gz'
        is_fastq = fastqutils.is_fastq(input_file)
        self.assertTrue(is_fastq)

    def test_is_fastq_with_fasta(self) -> None:
        """
        Tests the function that checks whether the input file is a FASTQ file using a FASTA file as input.
        :return: None
        """
        input_file = TestFastqUtils.test_file_dir / 'toy.fasta'
        is_fastq = fastqutils.is_fastq(input_file)
        self.assertFalse(is_fastq)


if __name__ == '__main__':
    unittest.main()
