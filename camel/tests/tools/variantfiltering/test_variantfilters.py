import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.variantfiltering.depthfilter import DepthFilter
from camel.app.tools.variantfiltering.distancefilter import DistanceFilter
from camel.app.tools.variantfiltering.mappingqualityfilter import MappingQualityFilter
from camel.app.tools.variantfiltering.snpqualityfilter import SnpQualityFilter
from camel.app.tools.variantfiltering.zscorefilter import ZScoreFilter


class TestVariantFiltering(CamelTestSuite):
    """
    Tests the variant filters.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('variant_calling')

    PATH_VCF_GZ_UNFILTERED = test_file_dir / 'unfiltered_variants-myco.vcf.gz'
    PATH_BAM = test_file_dir / 'alignment.bam'
    PATH_BED = test_file_dir / 'h37Rv.bed'

    def test_depth_filter(self) -> None:
        """
        Tests the depth filter.
        :return: None
        """
        depth_filter = DepthFilter()
        depth_filter.add_input_files({'VCF_GZ': [ToolIOFile(TestVariantFiltering.PATH_VCF_GZ_UNFILTERED)]})
        depth_filter.update_parameters(min_depth=20, min_forward_depth=2, min_reverse_depth=2)
        depth_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in depth_filter.tool_outputs)
        self.assertGreater(depth_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty")

    def test_depth_filter_soft_mask(self) -> None:
        """
        Tests the depth filter with soft masking.
        :return: None
        """
        depth_filter = DepthFilter()
        depth_filter.add_input_files({'VCF_GZ': [ToolIOFile(TestVariantFiltering.PATH_VCF_GZ_UNFILTERED)]})
        depth_filter.update_parameters(min_depth=80, min_forward_depth=5, min_reverse_depth=5, soft_filter=True)
        depth_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in depth_filter.tool_outputs)
        self.assertGreater(depth_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty")

    def test_mapping_quality_filter(self) -> None:
        """
        Tests the mapping quality filter.
        :return: None
        """
        mq_filter = MappingQualityFilter()
        mq_filter.add_input_files({'VCF_GZ': [ToolIOFile(TestVariantFiltering.PATH_VCF_GZ_UNFILTERED)]})
        mq_filter.update_parameters(min_mapping_quality=25)
        mq_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in mq_filter.tool_outputs)
        self.assertGreater(mq_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty")

    def test_snp_quality_filter(self) -> None:
        """
        Tests the SNP quality filter.
        :return: None
        """
        sq_filter = SnpQualityFilter()
        sq_filter.add_input_files({'VCF_GZ': [ToolIOFile(TestVariantFiltering.PATH_VCF_GZ_UNFILTERED)]})
        sq_filter.update_parameters(min_snp_quality=50)
        sq_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in sq_filter.tool_outputs)
        self.assertGreater(sq_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty")

    def test_distance_filter(self) -> None:
        """
        Tests the distance filter.
        :return: None
        """
        distance_filter = DistanceFilter()
        distance_filter.add_input_files({'VCF_GZ': [ToolIOFile(TestVariantFiltering.PATH_VCF_GZ_UNFILTERED)]})
        distance_filter.update_parameters(min_distance=10, keep_best=True)
        distance_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in distance_filter.tool_outputs)
        self.assertGreater(distance_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty")

    def test_distance_filter_soft_mask(self) -> None:
        """
        Tests the distance filter.
        :return: None
        """
        distance_filter = DistanceFilter()
        distance_filter.add_input_files({'VCF_GZ': [ToolIOFile(TestVariantFiltering.PATH_VCF_GZ_UNFILTERED)]})
        distance_filter.update_parameters(min_distance=10, keep_best=False, soft_filter=True)
        distance_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in distance_filter.tool_outputs)
        self.assertGreater(distance_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty")

    def test_zscore_filter(self) -> None:
        """
        Tests the Z-score filter.
        :return: None
        """
        zscore_filter = ZScoreFilter()
        zscore_filter.add_input_files({
            'VCF_GZ': [ToolIOFile(TestVariantFiltering.PATH_VCF_GZ_UNFILTERED)],
            'BAM': [ToolIOFile(TestVariantFiltering.PATH_BAM)]
        })
        zscore_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in zscore_filter.tool_outputs)
        self.assertGreater(zscore_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty")

    def test_zscore_filter_soft_mask(self) -> None:
        """
        Tests the Z-score filter.
        :return: None
        """
        zscore_filter = ZScoreFilter()
        zscore_filter.add_input_files({
            'VCF_GZ': [ToolIOFile(TestVariantFiltering.PATH_VCF_GZ_UNFILTERED)],
            'BAM': [ToolIOFile(TestVariantFiltering.PATH_BAM)]
        })
        zscore_filter.update_parameters(soft_filter=True)
        zscore_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in zscore_filter.tool_outputs)
        self.assertGreater(zscore_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty")


if __name__ == '__main__':
    unittest.main()
