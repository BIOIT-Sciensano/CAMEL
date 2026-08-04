import unittest
from pathlib import Path

from snakemake.iocontainers import OutputFiles

from camel.app.core.snakemake.snakemakeutils import get_rule_dir


class TestGetRuleDir(unittest.TestCase):
    """
    Tests the get_rule_dir utility function.
    """

    def test_single_output(self) -> None:
        """
        Returns the correct directory when there is a single output.
        :return: None
        """
        self.assertEqual(Path('bacmet/report'), get_rule_dir(['bacmet/report/html.iob']))

    def test_multiple_outputs_same_dir(self) -> None:
        """
        Returns the correct directory when all outputs share the same directory.
        :return: None
        """
        output = [
            'bacmet/report/html.iob',
            'bacmet/report/informs.io',
        ]
        self.assertEqual(Path('bacmet/report'), get_rule_dir(output))

    def test_multiple_outputs_different_dirs_raises(self) -> None:
        """
        Raises a ValueError when outputs span more than one directory.
        :return: None
        """
        output = [
            'bacmet/report/html.iob',
            'bacmet/tool/informs.io',
        ]
        with self.assertRaises(ValueError):
            get_rule_dir(output)

    def test_error_message_contains_dirs(self) -> None:
        """
        Checks that the ValueError message lists the conflicting directories.
        :return: None
        """
        output = [
            'bacmet/report/html.iob',
            'bacmet/tool/informs.io',
        ]
        with self.assertRaises(ValueError) as ctx:
            get_rule_dir(output)
        self.assertIn('bacmet/report', str(ctx.exception))
        self.assertIn('bacmet/tool', str(ctx.exception))

    def test_empty_output_raises(self) -> None:
        """
        Raises a ValueError when the output list is empty.
        :return: None
        """
        with self.assertRaises(ValueError):
            get_rule_dir([])

    def test_named_output_single(self) -> None:
        """
        Returns the correct directory when output is a single-entry Snakemake OutputFiles (named output).
        :return: None
        """
        output = OutputFiles(fromdict={'FASTA': 'bacmet/tool/fasta.io'})
        self.assertEqual(Path('bacmet/tool'), get_rule_dir(output))

    def test_named_output_multiple_same_dir(self) -> None:
        """
        Returns the correct directory when output is a multi-entry OutputFiles with all keys in the same directory.
        :return: None
        """
        output = OutputFiles(fromdict={
            'FASTA': 'bacmet/tool/fasta.io',
            'INFORMS': 'bacmet/tool/informs.io',
            'TSV': 'bacmet/tool/tsv.io',
        })
        self.assertEqual(Path('bacmet/tool'), get_rule_dir(output))

    def test_named_output_different_dirs_raises(self) -> None:
        """
        Raises a ValueError when a named OutputFiles has keys spread across different directories.
        :return: None
        """
        output = OutputFiles(fromdict={
            'FASTA': 'bacmet/tool/fasta.io',
            'HTML': 'bacmet/report/html.iob',
        })
        with self.assertRaises(ValueError):
            get_rule_dir(output)

    def test_named_output_empty_raises(self) -> None:
        """
        Raises a ValueError when an empty OutputFiles is passed.
        :return: None
        """
        with self.assertRaises(ValueError):
            get_rule_dir(OutputFiles())


if __name__ == '__main__':
    unittest.main()
