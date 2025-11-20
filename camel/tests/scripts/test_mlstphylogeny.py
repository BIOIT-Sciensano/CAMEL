import itertools
import unittest
from pathlib import Path

import pandas as pd

from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.scripts.mlstphylogeny.mainmlstphylogeny import main


class TestMLSTPhylogeny(CamelTestSuite):
    """
    Tests for the MLST phylogeny tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('mlst_phylogeny')
    dir_dataset_large_blast = test_file_dir / 'dataset_large_blast'
    dir_dataset_large_kma = test_file_dir / 'dataset_large_kma'
    dir_dataset_small_blast = test_file_dir / 'dataset_small_blast'
    dir_dataset_small_mist = test_file_dir / 'dataset_small_rapid'
    dir_dataset_small_blast_html = test_file_dir / 'dataset_html_small_blast'
    dir_dataset_small_blast_novel_alleles = test_file_dir / 'dataset_small_blast-novel_alleles'

    @staticmethod
    def generate_tsv_input_arguments(path_in: Path) -> list[str]:
        """
        Returns a list with the TSV input arguments.
        :param path_in: Directory with input files
        :return: List of input arguments
        """
        arguments = []
        for path_tsv in path_in.iterdir():
            if not path_tsv.is_file():
                continue
            arguments.extend(['--input-tsv', str(path_tsv), path_tsv.name])
        return arguments

    def test_mlst_phylogeny_small_dataset(self) -> None:
        """
        Tests the MLST phylogeny tool with standard options and a small dataset.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        result = cliutils.invoke(main, [
            *TestMLSTPhylogeny.generate_tsv_input_arguments(TestMLSTPhylogeny.dir_dataset_small_blast),
            '--output-html', str(output_html),
            '--output-dir', str(output_html.parent),
            '--dir-working', str(self.running_dir),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_html.stat().st_size, 0)

    def test_mlst_phylogeny_small_dataset_extra_outputs(self) -> None:
        """
        Tests the MLST phylogeny tool with standard options and a small dataset and extra output files.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_tsv_alleles = self.running_dir / 'out' / 'allele_matrix.tsv'
        output_tsv_dist = self.running_dir / 'out' / 'dist_matrix.tsv'
        output_nwk = self.running_dir / 'out' / 'tree.nwk'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        result = cliutils.invoke(main, [
            *TestMLSTPhylogeny.generate_tsv_input_arguments(TestMLSTPhylogeny.dir_dataset_small_blast),
            '--output-html', str(output_html),
            '--output-dir', str(output_html.parent),
            '--output-allele-matrix', str(output_tsv_alleles),
            '--output-dist-matrix', str(output_tsv_dist),
            '--output-tree', str(output_nwk),
            '--dir-working', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        for file in [output_html, output_tsv_alleles, output_tsv_dist, output_tsv_alleles]:
            self.assertGreater(file.stat().st_size, 0)

    def test_mlst_phylogeny_small_dataset_html(self) -> None:
        """
        Tests the MLST phylogeny tool with standard options and a small dataset and HTML input.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        input_tsv = list(itertools.chain.from_iterable(
            [['--input-html', str(dir_html)] for dir_html in TestMLSTPhylogeny.dir_dataset_small_blast_html.iterdir()]))
        args = input_tsv + [
            '--output-html', str(output_html),
            '--output-dir', str(output_html.parent),
            '--dir-working', str(self.running_dir),
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_html.stat().st_size, 0)

    def test_mlst_phylogeny_large_dataset(self) -> None:
        """
        Tests the MLST phylogeny tool with standard options and a large dataset.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        result = cliutils.invoke(
            main,
            [
                *TestMLSTPhylogeny.generate_tsv_input_arguments(
                    TestMLSTPhylogeny.dir_dataset_large_blast
                ),
                '--output-html', str(output_html),
                '--output-dir', str(output_html.parent),
                '--dir-working', str(self.running_dir),
            ],
        )
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_html.stat().st_size, 0)

    def test_mlst_phylogeny_large_dataset_kma(self) -> None:
        """
        Tests the MLST phylogeny tool with standard options and a large dataset and KMA-based detection.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        result = cliutils.invoke(
            main,
            [
                *TestMLSTPhylogeny.generate_tsv_input_arguments(
                    TestMLSTPhylogeny.dir_dataset_large_kma
                ),
                '--output-html', str(output_html),
                '--output-dir', str(output_html.parent),
                '--dir-working', str(self.running_dir),
                '--detection-method', 'kma',
            ],
        )
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_html.stat().st_size, 0)

    def test_mlst_phylogeny_identical_profiles(self) -> None:
        """
        Tests the MLST phylogeny tool with standard options and identical profiles for all datasets.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        input_tsv = list(sorted(TestMLSTPhylogeny.dir_dataset_small_blast.iterdir()))[0]
        args = [
            '--input-tsv', str(input_tsv), 'dataset_1.tsv',
            '--input-tsv', str(input_tsv), 'dataset_2.tsv',
            '--input-tsv', str(input_tsv), 'dataset_3.tsv',
            '--input-tsv', str(input_tsv), 'dataset_4.tsv',
            '--output-html', str(output_html),
            '--output-dir', str(output_html.parent),
            '--dir-working', str(self.running_dir),
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_html.stat().st_size, 0)

    def test_mlst_phylogeny_all_loci(self) -> None:
        """
        Tests the MLST phylogeny tool with standard options and a large dataset while keeping all loci.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        output_tsv_alleles = self.running_dir / 'out' / 'allele_matrix.tsv'
        result = cliutils.invoke(
            main,
            [
                *TestMLSTPhylogeny.generate_tsv_input_arguments(
                    TestMLSTPhylogeny.dir_dataset_large_blast
                ),
                '--output-html', str(output_html),
                '--output-dir', str(output_html.parent),
                '--output-allele-matrix', str(output_tsv_alleles),
                '--dir-working', str(self.running_dir),
                '--keep-all-loci',
            ],
        )
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_html.stat().st_size, 0)

        # Check if all loci are retained
        nb_loci_in = len(pd.read_table(next(iter(TestMLSTPhylogeny.dir_dataset_large_blast.iterdir()))).index)
        nb_loci_out = len(pd.read_table(output_tsv_alleles).columns) - 1
        self.assertEqual(nb_loci_in, nb_loci_out)

    def test_mlst_phylogeny_novel_alleles(self) -> None:
        """
        Tests the MLST phylogeny tool with novel alleles.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        result = cliutils.invoke(
            main,
            [
                *TestMLSTPhylogeny.generate_tsv_input_arguments(
                    TestMLSTPhylogeny.dir_dataset_small_blast_novel_alleles
                ),
                '--output-html', str(output_html),
                '--output-dir', str(output_html.parent),
                '--dir-working', str(self.running_dir),
                '--min-perc-samples', '50',
                '--min-perc-loci', '50',
            ],
        )
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_html.stat().st_size, 0)

    def test_mlst_phylogeny_novel_alleles_no_temp(self) -> None:
        """
        Tests the MLST phylogeny tool with novel alleles and no temporary alleles.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        result = cliutils.invoke(
            main,
            [
                *TestMLSTPhylogeny.generate_tsv_input_arguments(
                    TestMLSTPhylogeny.dir_dataset_small_blast_novel_alleles
                ),
                '--output-html', str(output_html),
                '--output-dir', str(output_html.parent),
                '--dir-working', str(self.running_dir),
                '--min-perc-samples', '50',
                '--min-perc-loci', '50',
                '--no-temp-allele-ids',
            ],
        )
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_html.stat().st_size, 0)

    def test_mlst_phylogeny_small_mist(self) -> None:
        """
        Tests the MLST phylogeny tool with standard options and a small dataset called with the 'mist' method.
        :return: None
        """
        output_html = self.running_dir / 'out' / 'report.html'
        output_html.parent.mkdir(exist_ok=True, parents=True)
        result = cliutils.invoke(
            main,
            [
                *TestMLSTPhylogeny.generate_tsv_input_arguments(
                    TestMLSTPhylogeny.dir_dataset_small_mist
                ),
                '--output-html', str(output_html),
                '--output-dir', str(output_html.parent),
                '--dir-working', str(self.running_dir),
                '--detection-method', 'mist',
                '--min-perc-loci', '75',
                '--min-perc-samples', '0',
            ],
        )
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_html.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
