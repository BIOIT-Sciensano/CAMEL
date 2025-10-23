import filecmp
from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import fileutils
from camel.app.core.cameltestsuite import CamelTestSuite


class TestFileUtils(CamelTestSuite):
    """
    Tests the Fastq utils module.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('components')

    def test_hash_file(self) -> None:
        """
        Tests the hash file function.
        :return: None
        """
        self.assertEqual(
            '258cf917517a15d335c42e9da9efdd6f157c88a22ab4a9aa94743dafe3c82fc1',
            fileutils.hash_file(TestFileUtils.test_file_dir / 'fq-fwd.fq'))

    def test_get_all_files(self) -> None:
        """
        Tests the get all files function.
        :return: None
        """
        self.assertCountEqual(
            fileutils.get_all_files(TestFileUtils.test_file_dir / 'hash_test'),
            [
                Path('/testdata/camel/components/hash_test/fq-fwd.fq'),
                Path('/testdata/camel/components/hash_test/fq-pe.fq'),
                Path('/testdata/camel/components/hash_test/fq-rev.fq'),
                Path('/testdata/camel/components/hash_test/fq-se.fq')
            ])

    def test_hash_directory(self) -> None:
        """
        Tests the hash directory function.
        :return: None
        """
        self.assertEqual(
            'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855',
            fileutils.hash_directory(TestFileUtils.test_file_dir / 'fq-fwd.fq')
        )

    def test_hash_value(self) -> None:
        """
        Tests the hash value function.
        :return: None
        """
        self.assertEqual(
            '332f61e385645bc5ad755f2ff6ccf8c4693f569a79f33427065731e9ff41b812',
            fileutils.hash_value('this_value_will_be_hashed')
        )

    def test_concatenate_files(self) -> None:
        """
        Tests the concatenate files function.
        :return: None
        """
        fileutils.concatenate_files(self.running_dir / 'concat.fq', [
            TestFileUtils.test_file_dir / 'fq-fwd.fq',
            TestFileUtils.test_file_dir / 'fq-rev.fq'])
        self.assertTrue(filecmp.cmp(
            TestFileUtils.test_file_dir / 'fq-fwd-rev-concatenated.fq',
            self.running_dir / 'concat.fq', shallow=False)
        )

    def test_concatenate_gzipped_files(self) -> None:
        """
        Tests the concatenate files function with gzipped files.
        :return: None
        """
        fileutils.concatenate_files(self.running_dir / 'concat.fq.gz', [
            TestFileUtils.test_file_dir / 'fq-fwd.fq.gz',
            TestFileUtils.test_file_dir / 'fq-rev.fq.gz'])
        cmd = Command(f'gunzip -c {self.running_dir / "concat.fq.gz"} > {self.running_dir / "unzipped.fq"}')
        cmd.run(self.running_dir)
        self.assertTrue(filecmp.cmp(TestFileUtils.test_file_dir / 'fq-fwd-rev-concatenated.fq', self.running_dir / 'unzipped.fq', shallow=False))
