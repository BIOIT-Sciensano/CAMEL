import unittest


from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.scripts.hybridassemblypipeline.mainhybridassemblypipeline import MainHybridAssemblyPipeline
from camel.tests import longRunningTest


class TestHybridAssemblyPipeline(CamelTestSuite):
    """
    Class to test the hybrid assembly pipeline.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('hybrid_assembly')
    FASTQ_1 = test_file_dir / 'ncbi_region_1.fastq'
    FASTQ_2 = test_file_dir / 'ncbi_region_2.fastq'
    FASTQ_SE = test_file_dir / 'ncbi_region_ont.fastq.gz'

    def test_sample_name(self) -> None:
        """
        Tests if the hybrid assembly pipeline extracts the sample name correctly.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'output.html'
        args = [
            '--output-html', str(path_report_out),
            '--fastq-pe', str(TestHybridAssemblyPipeline.FASTQ_1), str(TestHybridAssemblyPipeline.FASTQ_2),
            '--fastq-pe-names', 'illumina_1.fastq.gz', 'illumina_2.fastq.gz',
            '--fastq-se', str(TestHybridAssemblyPipeline.FASTQ_SE),
            '--fastq-se-name', 'ont.fastq.gz',
            '--working-dir', str(self.running_dir),
            '--expected-genome-size', '4.5m',
            '--ont-qual', 'nano-corr',
            '--input-type', 'hybrid'
        ]
        main = MainHybridAssemblyPipeline(args)
        self.assertEqual(main._sample_name, 'illumina')

    @longRunningTest()
    def test_hybrid_assembly(self) -> None:
        """
        Tests the hybrid assembly pipeline.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'output.html'
        args = [
            '--output-html', str(path_report_out),
            '--fastq-pe', str(TestHybridAssemblyPipeline.FASTQ_1), str(TestHybridAssemblyPipeline.FASTQ_2),
            '--fastq-se', str(TestHybridAssemblyPipeline.FASTQ_SE),
            '--working-dir', str(self.running_dir),
            '--input-type', 'hybrid',
            '--threads', '16',
            '--expected-genome-size', '4.5m',
            '--ont-qual', 'nano-corr',
            '--ont-basecalling-model', 'r1041_e82_400bps_sup_v5.0.0'
        ]
        main = MainHybridAssemblyPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_hybrid_assembly_with_unicycler(self) -> None:
        """
        Tests the hybrid assembly pipeline with Unicycler included.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'output.html'
        args = [
            '--output-html', str(path_report_out),
            '--fastq-pe', str(TestHybridAssemblyPipeline.FASTQ_1), str(TestHybridAssemblyPipeline.FASTQ_2),
            '--fastq-se', str(TestHybridAssemblyPipeline.FASTQ_SE),
            '--input-type', 'hybrid',
            '--working-dir', str(self.running_dir),
            '--threads', '16',
            '--expected-genome-size', '4.5m',
            '--ont-qual', 'nano-corr',
            '--unicycler'
        ]
        main = MainHybridAssemblyPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
