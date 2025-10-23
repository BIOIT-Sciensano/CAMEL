import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.scripts.pipelinecombine.mainpipelinecombine import MainPipelineCombine


class TestPipelineCombine(CamelTestSuite):
    """
    Tests the PipelineCombine tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('pipeline_combine')

    def test_pipeline_combine_staphylococcus(self) -> None:
        """
        Tests the pipeline combine main script.
        :return: None
        """
        path_tsv_out = self.running_dir / 'summary_combined.tsv'
        pipe_combine = MainPipelineCombine([
            *[str(x) for x in TestPipelineCombine.test_file_dir.glob('pipe_staphylococcus_*')],
            '--output', str(path_tsv_out)
        ])
        pipe_combine.run()
        self.assertGreater(path_tsv_out.stat().st_size, 0)

    def test_pipeline_combine_staphylococcus_with_formatted_genes(self) -> None:
        """
        Tests the pipeline combine main script.
        :return: None
        """
        path_tsv_out = self.running_dir / 'summary_combined.tsv'
        pipe_combine = MainPipelineCombine([
            *[str(x) for x in TestPipelineCombine.test_file_dir.glob('pipe_staphylococcus_*')],
            '--output', str(path_tsv_out),
            '--exclude', 'cgmlst-*,qc*,input_files,downsampling*,assembly*',
            '--gene-format', 'simple'
        ])
        pipe_combine.run()
        self.assertGreater(path_tsv_out.stat().st_size, 0)

    def test_pipeline_combine_staphylococcus_exclude(self) -> None:
        """
        Tests the pipeline combine main script.
        :return: None
        """
        path_tsv_out = self.running_dir / 'summary_combined.tsv'
        pipe_combine = MainPipelineCombine([
            *[str(x) for x in TestPipelineCombine.test_file_dir.glob('pipe_staphylococcus_*')],
            '--output', str(path_tsv_out),
            '--exclude', 'mlst-*,cgmlst-*,qc*,input_files,downsampling*,assembly*,vfdb*,trimming*',
            '--gene-format', 'locus_with_id'
        ])
        pipe_combine.run()
        self.assertGreater(path_tsv_out.stat().st_size, 0)

    def test_pipeline_combine_staphylococcus_include(self) -> None:
        """
        Tests the pipeline combine main script.
        :return: None
        """
        path_tsv_out = self.running_dir / 'summary_combined.tsv'
        pipe_combine = MainPipelineCombine([
            *[str(x) for x in TestPipelineCombine.test_file_dir.glob('pipe_staphylococcus_*')],
            '--output', str(path_tsv_out),
            '--include', 'sample,mlst-ST,spa_type,hits_ncbi_amr',
            '--gene-format', 'simple'
        ])
        pipe_combine.run()
        self.assertGreater(path_tsv_out.stat().st_size, 0)

    def test_pipeline_combine_staphylococcus_lre_finder(self) -> None:
        """
        Tests the pipeline combine main script.
        :return: None
        """
        path_tsv_out = self.running_dir / 'summary_combined.tsv'
        pipe_combine = MainPipelineCombine([
            *[str(x) for x in TestPipelineCombine.test_file_dir.glob('pipe_staphylococcus_*')],
            '--output', str(path_tsv_out),
            '--include', 'sample,mlst-ST,spa_type,lrefinder_genes*,hits_ncbi_amr',
            '--gene-format', 'simple'
        ])
        pipe_combine.run()
        self.assertGreater(path_tsv_out.stat().st_size, 0)

    def test_pipeline_combine_staphylococcus_group_by_gene(self) -> None:
        """
        Tests the pipeline combine main script with grouping by gene name.
        :return: None
        """
        path_tsv_out = self.running_dir / 'summary_combined.tsv'
        pipe_combine = MainPipelineCombine([
            *[str(x) for x in TestPipelineCombine.test_file_dir.glob('pipe_staphylococcus_*')],
            '--output', str(path_tsv_out),
            '--exclude', 'mlst-*,cgmlst-*,qc*,input_files,downsampling*,assembly*,vfdb*,trimming*',
            '--gene-format', 'locus_with_id',
            '--group-genes', 'gene'
        ])
        pipe_combine.run()
        self.assertGreater(path_tsv_out.stat().st_size, 0)

    def test_pipeline_combine_staphylococcus_group_by_allele(self) -> None:
        """
        Tests the pipeline combine main script with grouping by allele.
        :return: None
        """
        path_tsv_out = self.running_dir / 'summary_combined.tsv'
        pipe_combine = MainPipelineCombine([
            *[str(x) for x in TestPipelineCombine.test_file_dir.glob('pipe_staphylococcus_*')],
            '--output', str(path_tsv_out),
            '--exclude', 'mlst-*,cgmlst-*,qc*,input_files,downsampling*,assembly*,vfdb*,trimming*',
            '--gene-format', 'locus_with_id',
            '--group-genes', 'allele'
        ])
        pipe_combine.run()
        self.assertGreater(path_tsv_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
