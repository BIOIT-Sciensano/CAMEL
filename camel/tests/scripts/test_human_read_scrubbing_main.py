import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.utils import fastqutils
from camel.scripts.ncbihumanreadscrubber.mainncbihumanreadscrubber import MainNcbiHumanReadScrubber


class TestNcbiHumanReadScrubber(CamelTestSuite):
    """
    Tests the HRRT tool.
    _nh: files without human reads
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('ncbi_human_read_scrubber')

    path_fq_se = test_file_dir / 'Myco-DRR041783-ds_1_subset.fastq'
    path_fq_ont = test_file_dir / 'minion_reads.fastq'
    path_fq_ont_nh = test_file_dir / 'minion_no_human.fastq'
    files_list = [path_fq_se, path_fq_ont, path_fq_ont_nh]
    output = [True, True, False]

    datasets_pe = [
        {
            'fq': ['nh_reads_illumina_1.fastq', 'nh_reads_illumina_2.fastq'],
            'contains_human': False
        },
        {
            'fq': ['nh_reads_illumina_1.fastq.gz', 'nh_reads_illumina_2.fastq.gz'],
            'contains_human': False
        },
        {
            'fq': ['reads_illumina_1.fastq', 'reads_illumina_2.fastq'],
            'contains_human': True
        },
        {
            'fq': ['reads_illumina_1.fastq.gz', 'reads_illumina_2.fastq.gz'],
            'contains_human': True
        }
    ]

    datasets_ont = [
        {
            'fq': 'minion_reads.fastq',
            'contains_human': True
        },
        {
            'fq': 'minion_no_human.fastq',
            'contains_human': False
        },
    ]

    def test_scrubbing_illumina(self) -> None:
        """
        Tests the NCBI human read scrubbing pipeline on PE Illumina input.
        :return: None
        """
        for idx, dataset in enumerate(TestNcbiHumanReadScrubber.datasets_pe):
            # Create a novel working dir
            dir_working = self.running_dir / f'dataset_{idx}'
            dir_working.mkdir(parents=True, exist_ok=True)

            # Setup the output file
            path_report = dir_working / 'out' / 'report.html'

            # Run the pipeline
            args = (
                '--fastq-pe', *(str(TestNcbiHumanReadScrubber.test_file_dir / fq) for fq in dataset['fq']),
                '--output-html', str(path_report),
                '--output-dir', str(path_report.parent),
                '--working-dir', str(dir_working),
                '--output-tsv', "None",
                '--input-type', 'illumina',
                '--threads', '2'
            )
            main = MainNcbiHumanReadScrubber(args)
            main.run()

            # Check the output files
            self.assertGreater(path_report.stat().st_size, 0)
            for fq_idx, ext in enumerate(['_R1.fastq.gz', '_R2.fastq.gz']):
                path_fq_scrubbed = next(path_report.parent.glob(f"*{ext}"))
                self.assertGreater(path_fq_scrubbed.stat().st_size, 0)

                if dataset['contains_human'] is False:
                    # Check whether all reads were retained (no human reads)
                    self.assertEqual(
                        fastqutils.count_reads(TestNcbiHumanReadScrubber.test_file_dir / dataset['fq'][fq_idx]),
                        fastqutils.count_reads(path_fq_scrubbed),
                        "Input and output file should contain the same number of reads."
                    )
                else:
                    # Check whether reads were removed (human reads)
                    self.assertGreater(
                        fastqutils.count_reads(TestNcbiHumanReadScrubber.test_file_dir / dataset['fq'][fq_idx]),
                        fastqutils.count_reads(path_fq_scrubbed),
                        "No reads removed from dataset containing human reads."
                    )

    def test_scrubbing_illumina_with_export(self) -> None:
        """
        Tests the NCBI human read scrubbing pipeline on PE Illumina input.
        The removed reads are exported.
        :return: None
        """
        for idx, dataset in enumerate(TestNcbiHumanReadScrubber.datasets_pe):
            # Create a novel working dir
            dir_working = self.running_dir / f'dataset_{idx}'
            dir_working.mkdir(parents=True, exist_ok=True)

            # Setup the output file
            path_report = dir_working / 'out' / 'report.html'

            # Run the pipeline
            args = (
                '--fastq-pe', *(str(TestNcbiHumanReadScrubber.test_file_dir / fq) for fq in dataset['fq']),
                '--output-html', str(path_report),
                '--output-dir', str(path_report.parent),
                '--working-dir', str(dir_working),
                '--output-tsv', "None",
                '--input-type', 'illumina',
                '--threads', '2',
                '--export-removed-reads'
            )
            main = MainNcbiHumanReadScrubber(args)
            main.run()

            # Check the output files
            for fq_idx in (1, 2):
                if dataset['contains_human'] is False:
                    with self.assertRaises(StopIteration):
                        next(path_report.parent.glob(f"*-removed_R{fq_idx}.fastq.gz"))
                else:
                    try:
                        path_fq = next(path_report.parent.glob(f"*-removed_R{fq_idx}.fastq.gz"))
                    except StopIteration:
                        paths_fq = (fq.name for fq in path_report.parent.glob('*.fastq.gz'))
                        raise FileNotFoundError(f"File with removed reads _R{fq_idx} not found ({', '.join(paths_fq)})")
                    self.assertGreater(fastqutils.count_reads(path_fq), 0)

    def test_scrubbing_ont(self) -> None:
        """
        Tests the NCBI human read scrubbing pipeline on SE ONT input.
        :return: None
        """
        for idx, dataset in enumerate(TestNcbiHumanReadScrubber.datasets_ont):
            # Create a novel working dir
            dir_working = self.running_dir / f'dataset_{idx}'
            dir_working.mkdir(parents=True, exist_ok=True)

            # Setup the output file
            path_report = dir_working / 'out' / 'report.html'

            # Run the pipeline
            args = (
                '--fastq-se', str(TestNcbiHumanReadScrubber.test_file_dir / dataset['fq']),
                '--output-html', str(path_report),
                '--output-dir', str(path_report.parent),
                '--working-dir', str(dir_working),
                '--output-tsv', "None",
                '--input-type', 'ont',
                '--threads', '2'
            )
            main = MainNcbiHumanReadScrubber(args)
            main.run()

            # Check the output files
            self.assertGreater(path_report.stat().st_size, 0)
            path_fq_scrubbed = next(path_report.parent.glob("*-scrubbed.fastq.gz"))
            self.assertGreater(path_fq_scrubbed.stat().st_size, 0)

            if dataset['contains_human'] is False:
                # Check whether all reads were retained (no human reads)
                self.assertEqual(
                    fastqutils.count_reads(TestNcbiHumanReadScrubber.test_file_dir / dataset['fq']),
                    fastqutils.count_reads(path_fq_scrubbed),
                    "Input and output file should contain the same number of reads."
                )
            else:
                # Check whether reads were removed (human reads)
                self.assertGreater(
                    fastqutils.count_reads(TestNcbiHumanReadScrubber.test_file_dir / dataset['fq']),
                    fastqutils.count_reads(path_fq_scrubbed),
                    "No reads removed from dataset containing human reads."
                )

    def test_scrubbing_ont_with_export(self) -> None:
        """
        Tests the NCBI human read scrubbing pipeline on SE ONT input.
        The removed reads are exported.
        :return: None
        """
        for idx, dataset in enumerate(TestNcbiHumanReadScrubber.datasets_ont):
            # Create a novel working dir
            dir_working = self.running_dir / f'dataset_{idx}'
            dir_working.mkdir(parents=True, exist_ok=True)

            # Setup the output file
            path_report = dir_working / 'out' / 'report.html'

            # Run the pipeline
            args = (
                '--fastq-se', str(TestNcbiHumanReadScrubber.test_file_dir / dataset['fq']),
                '--output-html', str(path_report),
                '--output-dir', str(path_report.parent),
                '--working-dir', str(dir_working),
                '--output-tsv', "None",
                '--input-type', 'ont',
                '--threads', '2',
                '--export-removed-reads'
            )
            main = MainNcbiHumanReadScrubber(args)
            main.run()

            # Check the output files
            self.assertGreater(path_report.stat().st_size, 0)
            path_fq_scrubbed = next(path_report.parent.glob("*-scrubbed.fastq.gz"))

            if dataset['contains_human'] is False:
                with self.assertRaises(StopIteration):
                    next(path_report.parent.glob("*-removed.fastq.gz"))
            else:
                path_removed = next(path_report.parent.glob("*-removed.fastq.gz"))
                # Check whether reads were removed (human reads)
                self.assertEqual(
                    fastqutils.count_reads(TestNcbiHumanReadScrubber.test_file_dir / dataset['fq']),
                    fastqutils.count_reads(path_fq_scrubbed) + fastqutils.count_reads(path_removed),
                    "Number of reads should match."
                )

        # def test_scrubbing_paired(self) -> None:
    #     """
    #     Tests the NCBI human read scrubbing standalone pipeline with illumina files.
    #     :return: None
    #     """
    #     path_report_out = self.running_dir / 'out' / 'report.html'
    #     path_removed_reads = self.running_dir / 'out' / 'human_read_scrubbing'
    #     output_file_1 = self.running_dir / 'human_read_scrubbing' / 'fastq_pe' / 'output' / 'fastq_removed_1.fastq.gz'
    #     output_file_2 = self.running_dir / 'human_read_scrubbing' / 'fastq_pe' / 'output' / 'fastq_removed_2.fastq.gz'
    #     # _nh: files without human reads
    #     hum_reads = ['', 'nh_', '', 'nh_']
    #     extension = ['', '', '.gz', '.gz']
    #
    #     for hr, ext in zip(hum_reads, extension):
    #         args = [
    #             '--fastq-pe',
    #             str(TestNcbiHumanReadScrubber.test_file_dir / f'{hr}reads_illumina_1.fastq{ext}'),
    #             str(TestNcbiHumanReadScrubber.test_file_dir / f'{hr}reads_illumina_2.fastq{ext}'),
    #             '--output-html', str(path_report_out),
    #             '--output-dir', str(path_report_out.parent),
    #             '--working-dir', str(self.running_dir),
    #             '--output-tsv', "None",
    #             '--input-type', 'illumina',
    #             '--threads', '2',
    #             '--export-removed-reads'
    #         ]
    #         main = MainNcbiHumanReadScrubber(args)
    #         main.run()
    #         self.assertGreater(path_report_out.stat().st_size, 0)
    #         self.assertGreater(output_file_1.stat().st_size, 0) and self.assertGreater(output_file_2.stat().st_size, 0) if hr == '' else not self.assertTrue(path_removed_reads.exists())
    #
    # def test_scrubbing_paired_wo_removed_reads(self) -> None:
    #     """
    #     Tests the NCBI human read scrubbing standalone pipeline with illumina files.
    #     :return: None
    #     """
    #     path_report_out = self.running_dir / 'out' / 'report.html'
    #     # _nh: files without human reads
    #     hum_reads = ['', 'nh_', '', 'nh_']
    #     extension = ['', '', '.gz', '.gz']
    #
    #     for hr, ext in zip(hum_reads, extension):
    #         args = [
    #             '--fastq-pe',
    #             str(TestNcbiHumanReadScrubber.test_file_dir / f'{hr}reads_illumina_1.fastq{ext}'),
    #             str(TestNcbiHumanReadScrubber.test_file_dir / f'{hr}reads_illumina_2.fastq{ext}'),
    #             '--output-html', str(path_report_out),
    #             '--output-dir', str(path_report_out.parent),
    #             '--working-dir', str(self.running_dir),
    #             '--output-tsv', "None",
    #             '--input-type', 'illumina',
    #             '--threads', '2',
    #         ]
    #         main = MainNcbiHumanReadScrubber(args)
    #         main.run()
    #         self.assertGreater(path_report_out.stat().st_size, 0)
    #
    # def test_scrubbing_fasta(self) -> None:
    #     """
    #     Tests the NCBI human read scrubbing standalone pipeline with fasta files.
    #     :return: None
    #     """
    #     path_report_html = self.running_dir / 'out' / 'report.html'
    #     path_removed_reads = self.running_dir / 'out' / 'human_read_scrubbing'
    #     human_reads = ['', 'no_']
    #     output_file = self.running_dir / 'human_read_scrubbing' / 'fasta' / 'output' / 'fastq_removed.fasta'
    #
    #     for hum in human_reads:
    #         args = [
    #             '--fasta',
    #             str(TestNcbiHumanReadScrubber.test_file_dir / f'bacteria_seq_{hum}hum.fasta'),
    #             '--output-html', str(path_report_html),
    #             '--output-dir', str(path_report_html.parent),
    #             '--working-dir', str(self.running_dir),
    #             '--output-tsv', "None",
    #             '--input-type', 'fasta',
    #             '--threads', '2',
    #             '--export-removed-reads'
    #         ]
    #         main = MainNcbiHumanReadScrubber(args)
    #         main.run()
    #         self.assertGreater(path_report_html.stat().st_size, 0)
    #         self.assertGreater(output_file.stat().st_size, 0) if hum == '' else not self.assertTrue(path_removed_reads.exists())


if __name__ == '__main__':
    unittest.main()
