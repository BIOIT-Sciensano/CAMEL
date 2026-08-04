import unittest
from pathlib import Path

import pysam
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import vcfutils

from camel.app.core.cameltestsuite import CamelTestSuite
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

    @staticmethod
    def create_mock_vcf(
        path: Path, variants: list[tuple[str, int, float, str]]
    ) -> Path:
        """
        Creates a minimal indexed VCF.gz file for testing.
        :param path: Output .vcf.gz path
        :param variants: List of tuples: (chrom, pos, qual, filter_value)
        :returns: Path to the created VCF.gz
        """
        header = pysam.VariantHeader()

        header.add_meta("fileformat", value="VCFv4.2")
        header.contigs.add("chr1", length=1000)
        header.contigs.add("chr2", length=500)

        header.add_sample("SAMPLE")

        with pysam.VariantFile(path, "wz", header=header) as vcf:
            for chrom, pos, qual, filt in variants:
                record = vcf.new_record(
                    contig=chrom,
                    start=pos - 1,  # pysam uses 0-based starts
                    stop=pos,
                    alleles=("A", "T"),
                    qual=qual,
                    filter=filt,
                )
                vcf.write(record)
        pysam.tabix_index(str(path), preset="vcf", force=True)
        return path

    def test_depth_filter(self) -> None:
        """
        Tests the depth filter.
        :return: None
        """
        depth_filter = DepthFilter()
        depth_filter.add_input_files(
            {'VCF_GZ': [ToolIOFile(TestVariantFiltering.PATH_VCF_GZ_UNFILTERED)]}
        )
        depth_filter.update_parameters(
            min_depth=20, min_forward_depth=2, min_reverse_depth=2
        )
        depth_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in depth_filter.tool_outputs)
        self.assertGreater(
            depth_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty"
        )

    def test_depth_filter_soft_mask(self) -> None:
        """
        Tests the depth filter with soft masking.
        :return: None
        """
        depth_filter = DepthFilter()
        depth_filter.add_input_files(
            {'VCF_GZ': [ToolIOFile(TestVariantFiltering.PATH_VCF_GZ_UNFILTERED)]}
        )
        depth_filter.update_parameters(
            min_depth=80, min_forward_depth=5, min_reverse_depth=5, soft_filter=True
        )
        depth_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in depth_filter.tool_outputs)
        self.assertGreater(
            depth_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty"
        )

    def test_mapping_quality_filter(self) -> None:
        """
        Tests the mapping quality filter.
        :return: None
        """
        mq_filter = MappingQualityFilter()
        mq_filter.add_input_files(
            {'VCF_GZ': [ToolIOFile(TestVariantFiltering.PATH_VCF_GZ_UNFILTERED)]}
        )
        mq_filter.update_parameters(min_mapping_quality=25)
        mq_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in mq_filter.tool_outputs)
        self.assertGreater(
            mq_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty"
        )

    def test_snp_quality_filter(self) -> None:
        """
        Tests the SNP quality filter.
        :return: None
        """
        sq_filter = SnpQualityFilter()
        sq_filter.add_input_files(
            {'VCF_GZ': [ToolIOFile(TestVariantFiltering.PATH_VCF_GZ_UNFILTERED)]}
        )
        sq_filter.update_parameters(min_snp_quality=50)
        sq_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in sq_filter.tool_outputs)
        self.assertGreater(
            sq_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty"
        )

    def test_distance_filter(self) -> None:
        """
        Tests the distance filter.
        :return: None
        """
        distance_filter = DistanceFilter()
        distance_filter.add_input_files(
            {'VCF_GZ': [ToolIOFile(TestVariantFiltering.PATH_VCF_GZ_UNFILTERED)]}
        )
        distance_filter.update_parameters(min_distance=10, keep_best=True)
        distance_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in distance_filter.tool_outputs)
        self.assertGreater(
            distance_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty"
        )

    def test_distance_filter_soft_mask(self) -> None:
        """
        Tests the distance filter.
        :return: None
        """
        distance_filter = DistanceFilter()
        distance_filter.add_input_files(
            {'VCF_GZ': [ToolIOFile(TestVariantFiltering.PATH_VCF_GZ_UNFILTERED)]}
        )
        distance_filter.update_parameters(
            min_distance=10, keep_best=False, soft_filter=True
        )
        distance_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in distance_filter.tool_outputs)
        self.assertGreater(
            distance_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty"
        )

    def test_distance_filter_with_mock(self) -> None:
        """
        Tests the depth filter with a mock VCF file to test edge cases.
        :return: None
        """
        path_mock = self.running_dir / 'mock.vcf.gz'
        TestVariantFiltering.create_mock_vcf(
            path_mock,
            [
                ("chr1", 100, 40, 'PASS'),
                (
                    "chr1",
                    150,
                    10,
                    'PASS',
                ),  # Should be removed (less than 100 bp and lower quality)
                ("chr1", 250, 40, 'PASS'),
                ("chr2", 100, 40, 'PASS'),
                (
                    "chr2",
                    300,
                    20,
                    'PASS',
                ),  # Should be removed (less than 10 bp and lower quality)
                ("chr2", 320, 40, 'PASS'),
            ],
        )
        distance_filter = DistanceFilter()
        distance_filter.add_input_files({'VCF_GZ': [ToolIOFile(path_mock)]})
        distance_filter.update_parameters(min_distance=100, keep_best=True)
        distance_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in distance_filter.tool_outputs)
        self.assertGreater(
            distance_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty"
        )
        self.assertEqual(
            vcfutils.count_variants(distance_filter.tool_outputs['VCF_GZ'][0].path), 4
        )

    def test_zscore_filter(self) -> None:
        """
        Tests the Z-score filter.
        :return: None
        """
        zscore_filter = ZScoreFilter()
        zscore_filter.add_input_files(
            {
                'VCF_GZ': [ToolIOFile(TestVariantFiltering.PATH_VCF_GZ_UNFILTERED)],
                'BAM': [ToolIOFile(TestVariantFiltering.PATH_BAM)],
            }
        )
        zscore_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in zscore_filter.tool_outputs)
        self.assertGreater(
            zscore_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty"
        )

    def test_zscore_filter_soft_mask(self) -> None:
        """
        Tests the Z-score filter.
        :return: None
        """
        zscore_filter = ZScoreFilter()
        zscore_filter.add_input_files(
            {
                'VCF_GZ': [ToolIOFile(TestVariantFiltering.PATH_VCF_GZ_UNFILTERED)],
                'BAM': [ToolIOFile(TestVariantFiltering.PATH_BAM)],
            }
        )
        zscore_filter.update_parameters(soft_filter=True)
        zscore_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in zscore_filter.tool_outputs)
        self.assertGreater(
            zscore_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty"
        )


if __name__ == '__main__':
    unittest.main()
