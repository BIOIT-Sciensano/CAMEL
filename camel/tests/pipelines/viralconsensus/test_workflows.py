import logging
import unittest
from pathlib import Path

from camel.app.config import config
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.scriptutils.fastqinput import FastqInput
from camel.app.core.utils import vcfutils
from camel.app.wrappers.trimmingilluminawrapper import TrimmingIlluminaWrapper
from camel.scripts.viralconsensuspipeline.workflows.callvariants import CallVariants
from camel.scripts.viralconsensuspipeline.workflows.filtervariants import FilterVariants
from camel.scripts.viralconsensuspipeline.workflows.readmappingworkflow import ReadMappingWorkflow
from camel.scripts.viralconsensuspipeline.workflows.segmentdownsampling import SegmentDownsamplingWorkflow
from camel.tests import longRunningTest


class TestWorkflows(CamelTestSuite):
    """
    Tests for the sub-workflows for the viral consensus pipeline.
    """

    def test_read_mapping_workflow_illumina(self) -> None:
        """
        Tests the read mapping workflow.
        :return: None
        """
        test_file_dir_illumina = CamelTestSuite.get_test_file_dir('read_mapping')
        fasta_ref = test_file_dir_illumina / 'reference.fasta'
        fq_fwd = test_file_dir_illumina / 'reads_1.fastq'
        fq_rev = test_file_dir_illumina / 'reads_2.fastq'
        fastq_input = FastqInput('illumina', [ToolIOFile(fq_fwd), ToolIOFile(fq_rev)], is_pe=True)
        workflow = ReadMappingWorkflow(self.running_dir)
        output = workflow.run(fastq_input, fasta_ref)
        self.assertTrue(output.path_bam.exists())
        self.assertGreater(output.path_bam.stat().st_size, 0)
        self.assertGreater(len(output.informs), 0)

    def test_read_mapping_workflow_illumina_depth_cutoff(self) -> None:
        """
        Tests the read mapping workflow.
        :return: None
        """
        test_file_dir_illumina = CamelTestSuite.get_test_file_dir('read_mapping')
        fasta_ref = test_file_dir_illumina / 'reference.fasta'
        fq_fwd = test_file_dir_illumina / 'reads_1.fastq'
        fq_rev = test_file_dir_illumina / 'reads_2.fastq'
        fastq_input = FastqInput('illumina', [ToolIOFile(fq_fwd), ToolIOFile(fq_rev)], is_pe=True)
        workflow = ReadMappingWorkflow(self.running_dir)
        output = workflow.run(fastq_input, fasta_ref, gap_depth_cutoff=1, gap_len_cutoff=3)
        self.assertTrue(output.path_bam.exists())
        self.assertGreater(output.path_bam.stat().st_size, 0)
        self.assertGreater(len(output.informs), 0)

    def test_read_mapping_workflow_illumina_trimmed(self) -> None:
        """
        Tests the read mapping workflow on trimmed reads.
        :return: None
        """
        # Input data
        test_file_dir_illumina = CamelTestSuite.get_test_file_dir('read_mapping')
        fasta_ref = test_file_dir_illumina / 'reference.fasta'
        fq_fwd = test_file_dir_illumina / 'reads_1.fastq'
        fq_rev = test_file_dir_illumina / 'reads_2.fastq'

        # Trim reads
        dir_trimming = self.running_dir / 'trim'
        wrapper = TrimmingIlluminaWrapper(dir_trimming)
        wrapper.run([fq_fwd, fq_rev], threads=4)
        fastq_input = FastqInput(
            read_type='illumina',
            pe=wrapper.output.trimmed_reads_pe,
            se_fwd=wrapper.output.trimmed_reads_se_fwd,
            se_rev=wrapper.output.trimmed_reads_se_rev,
            is_trimmed=True, is_pe=True)

        # Run read mapping workflow
        workflow = ReadMappingWorkflow(self.running_dir)
        output = workflow.run(fastq_input, fasta_ref, gap_len_cutoff=3, gap_depth_cutoff=12)
        self.assertTrue(output.path_bam.exists())
        self.assertGreater(output.path_bam.stat().st_size, 0)
        self.assertGreater(len(output.informs), 0)

    def test_read_mapping_workflow_ont(self) -> None:
        """
        Tests the read mapping workflow.
        :return: None
        """
        test_file_dir_ont = CamelTestSuite.get_test_file_dir('minion')
        fq_ont = test_file_dir_ont / 'fastq_minion_stec.fastq'
        fasta_ref_ont = test_file_dir_ont / 'NC_002695.2.fasta'
        fastq_input = FastqInput('nanopore', se=[ToolIOFile(fq_ont)], is_pe=False)
        workflow = ReadMappingWorkflow(self.running_dir)
        output = workflow.run(fastq_input, fasta_ref_ont)
        self.assertTrue(output.path_bam.exists())
        self.assertGreater(output.path_bam.stat().st_size, 0)
        self.assertGreater(len(output.informs), 0)

    def test_read_mapping_workflow_ont_no_depth(self) -> None:
        """
        Tests the read mapping workflow.
        :return: None
        """
        test_file_dir_ont = CamelTestSuite.get_test_file_dir('minion')
        fq_ont = test_file_dir_ont / 'fastq_minion_stec.fastq'
        fasta_ref_ont = test_file_dir_ont / 'H1N1_HA.fasta'
        fastq_input = FastqInput('nanopore', se=[ToolIOFile(fq_ont)], is_pe=False)
        workflow = ReadMappingWorkflow(self.running_dir)
        output = workflow.run(fastq_input, fasta_ref_ont)
        self.assertTrue(output.path_bam.exists())
        self.assertGreater(output.path_bam.stat().st_size, 0)
        self.assertGreater(len(output.informs), 0)

    def test_call_variants_ont_bcftools(self) -> None:
        """
        Runs the variant calling workflow with ONT data using bcftools.
        :return: None
        """
        dir_test = CamelTestSuite.get_test_file_dir('clair3')
        workflow = CallVariants(self.running_dir)
        output = workflow.run(dir_test / 'bsubtilis_ont.bam', dir_test / 'bsubtilis.fa', 'nanopore', 'bcftools', {})
        self.assertTrue(output.path_vcf.exists())
        self.assertGreater(vcfutils.count_variants(output.path_vcf), 0)
        logging.info(f'Stats: {output.stats}')

    def test_call_variants_illumina_bcftools(self) -> None:
        """
        Runs the variant calling workflow with ONT data using bcftools.
        :return: None
        """
        dir_test = CamelTestSuite.get_test_file_dir('clair3')
        workflow = CallVariants(self.running_dir)
        output = workflow.run(
            dir_test / 'bsubtilis_illumina.bam', dir_test / 'bsubtilis.fa', 'illumina', 'bcftools', {})
        self.assertTrue(output.path_vcf.exists())
        self.assertGreater(vcfutils.count_variants(output.path_vcf), 0)
        self.assertGreater(len(output.informs), 0)
        logging.info(f'Stats: {output.stats}')

    def test_filter_variants_bcftools(self) -> None:
        """
        Tests the filtering of variants using bcftools.
        :return: None
        """
        dir_test = CamelTestSuite.get_test_file_dir('variant_calling')
        workflow = FilterVariants(self.running_dir)
        filters = {'min_af': 0.5, 'min_dp': 20, 'min_qual': 25}
        out = workflow.run(dir_test / 'unfiltered_variants-myco.vcf', 'bcftools', filters)
        self.assertTrue(out.path_vcf.exists())
        self.assertGreater(len(out.informs), 0)
        self.assertGreater(vcfutils.count_variants(out.path_vcf), 0)

    @longRunningTest()
    def test_call_variants_ont_clair3(self) -> None:
        """
        Runs the variant calling workflow with ONT data using bcftools.
        :return: None
        """
        dir_test = CamelTestSuite.get_test_file_dir('clair3')
        model = Path(config.dir_db, 'clair3', 'models', 'ont')
        workflow = CallVariants(self.running_dir)
        output = workflow.run(
            dir_test / 'bsubtilis_ont.bam', dir_test / 'bsubtilis.fa', 'nanopore', 'clair3', {'model': model})
        self.assertTrue(output.path_vcf.exists())
        self.assertGreater(vcfutils.count_variants(output.path_vcf), 0)
        self.assertGreater(len(output.informs), 0)
        logging.info(f'Stats: {output.stats}')

    @longRunningTest()
    def test_call_variants_illumina_clair3(self) -> None:
        """
        Runs the variant calling workflow with ONT data using bcftools.
        :return: None
        """
        dir_test = CamelTestSuite.get_test_file_dir('clair3')
        model = Path(config.dir_db, 'clair3', 'models', 'ilmn')
        workflow = CallVariants(self.running_dir)
        output = workflow.run(
            dir_test / 'bsubtilis_illumina.bam', dir_test / 'bsubtilis.fa', 'illumina', 'clair3', {'model': model})
        self.assertTrue(output.path_vcf.exists())
        self.assertGreater(vcfutils.count_variants(output.path_vcf), 0)
        self.assertGreater(len(output.informs), 0)
        logging.info(f'Stats: {output.stats}')

    def test_filter_variants_clair3(self) -> None:
        """
        Tests the filtering of variants using Clair3.
        :return: None
        """
        dir_test = CamelTestSuite.get_test_file_dir('variant_calling')
        workflow = FilterVariants(self.running_dir)
        filters = {'min_af': 0.5, 'min_dp': 20, 'min_qual': 25}
        out = workflow.run(dir_test / 'unfiltered_variants-myco-clair3.vcf.gz', 'clair3', filters)
        self.assertTrue(out.path_vcf.exists())
        self.assertGreater(len(out.informs), 0)
        self.assertGreater(vcfutils.count_variants(out.path_vcf), 0)

    def test_segment_downsampling_workflow(self) -> None:
        """
        Tests the segment downsampling workflow.
        """
        dir_test = CamelTestSuite.get_test_file_dir('pipelines', 'viral_consensus')
        path_json = dir_test / 'mapping_ilmn_infl_a.json'
        path_bam = dir_test / 'mapping_ilmn_infl_a.bam'
        workflow = SegmentDownsamplingWorkflow(self.running_dir)
        out = workflow.run(path_bam, 'illumina', path_json, 100, threads=8)
        self.assertTrue(all(io.path.exists() for io in out.fq_out.pe))
        self.assertGreater(len(out.informs), 0)

    def test_segment_downsampling_workflow_primer_clipping(self) -> None:
        """
        Tests the segment downsampling workflow with primer clipping.
        """
        # Run clipping
        dir_test = CamelTestSuite.get_test_file_dir('pipelines', 'viral_consensus')
        path_json = dir_test / 'mapping_ilmn_sars_cov_2.json'
        path_bam = dir_test / 'mapping_sars_cov_2.bam'
        path_bed = dir_test / 'mapping_sars_cov_2-primers.bed'
        workflow = SegmentDownsamplingWorkflow(self.running_dir)
        out = workflow.run(path_bam, 'illumina', path_json, 100, bed_primers=path_bed, threads=8)
        self.assertTrue(all(io.path.exists() for io in out.fq_out.pe))
        self.assertGreater(len(out.informs), 0)


if __name__ == '__main__':
    unittest.main()
