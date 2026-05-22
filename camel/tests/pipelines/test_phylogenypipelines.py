from pathlib import Path

from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.scripts.snpphylogeny.maincfsanphylo import main as main_cfsan
from camel.scripts.snpphylogeny.mainsamtoolsphylo import main as main_samtools
from camel.tests import longRunningTest


class TestSnpPhylogenyPipelines(CamelTestSuite):
    """
    Tests the SNP phylogeny pipelines.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('phylogeny', 'lm_subset')
    reference_fasta = test_file_dir / 'reference.fasta'

    @staticmethod
    def __get_samples(gzipped: bool = False) -> list[str]:
        """
        Returns the input samples.
        :return: Input samples
        """
        samples = []
        pattern = '*_1.fastq' + ('.gz' if gzipped else '')
        for index, fq_fwd in enumerate(TestSnpPhylogenyPipelines.test_file_dir.glob(pattern), 1):
            fq_rev = fq_fwd.parent / fq_fwd.name.replace('_1', '_2')
            samples.extend([
                '--sample', f'sample_{index}', fq_fwd.name, str(fq_fwd), fq_rev.name, str(fq_rev)])
        return samples

    @longRunningTest()
    def test_samtools_phylogeny(self) -> None:
        """
        Tests the Samtools Phylogeny pipeline.
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        output_file_fasta = Path(self.running_dir) / 'report' / 'snps.fasta'
        result = cliutils.invoke(main_samtools, TestSnpPhylogenyPipelines.__get_samples() + [
            '--reference', str(TestSnpPhylogenyPipelines.reference_fasta),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--output-fasta', str(output_file_fasta),
            '--working-dir', str(self.running_dir),
            '--trim-reads',
            '--include-ref',
            '--min-total-depth', '1',
            '--min-fwd-depth', '0',
            '--min-rev-depth', '0',
            '--min-distance', '2'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_file_report.stat().st_size, 0)
        self.assertGreater(output_file_fasta.stat().st_size, 0)

    @longRunningTest()
    def test_samtools_phylogeny_with_masking(self) -> None:
        """
        Tests the Samtools Phylogeny pipeline.
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        result = cliutils.invoke(main_samtools, TestSnpPhylogenyPipelines.__get_samples() + [
            '--reference', str(TestSnpPhylogenyPipelines.reference_fasta),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--trim-reads', '--include-ref',
            '--min-total-depth', '1',
            '--min-fwd-depth', '0',
            '--min-rev-depth', '0',
            '--min-distance', '2',
            '--soft-filter'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_file_report.stat().st_size, 0)

    @longRunningTest()
    def test_cfsan_phylogeny(self) -> None:
        """
        Tests the CFSAN Phylogeny pipeline.
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        output_file_fasta = Path(self.running_dir) / 'report' / 'snps.fasta'
        result = cliutils.invoke(main_cfsan, TestSnpPhylogenyPipelines.__get_samples() + [
            '--reference', str(TestSnpPhylogenyPipelines.reference_fasta),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--output-fasta', str(output_file_fasta),
            '--working-dir', str(self.running_dir),
            '--trim-reads',
            '--missing-data', 'complete_deletion',
            '--branch-swap', 'none',
            '--site-cov-cutoff', '50',
            '--bootstraps', '10'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_file_report.stat().st_size, 0)

    @longRunningTest()
    def test_cfsan_phylogeny_gz_no_trim(self) -> None:
        """
        Tests the CFSAN Phylogeny pipeline with gzipped input files and no trimming enabled.
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        result = cliutils.invoke(main_cfsan, TestSnpPhylogenyPipelines.__get_samples() + [
            '--reference', str(TestSnpPhylogenyPipelines.reference_fasta),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--missing-data', 'complete_deletion',
            '--ml-method', 'spr3',
            '--site-cov-cutoff', '50',
            '--bootstraps', '10'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_file_report.stat().st_size, 0)
