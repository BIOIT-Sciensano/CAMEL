import unittest

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.scriptutils.basepipe.fastqinput import FastqInput
from camel.app.wrappers.assemblywrapper import AssemblyWrapper


class TestWorkflowAssembly(CamelTestSuite):
    """
    Tests the Illumina trimming workflow.
    """
    # Input files (Gene detection)
    test_file_dir = CamelTestSuite.get_test_file_dir('workflows')

    # Input files
    fastq_pe = [
        test_file_dir / 'assembly' / 'reads_illumina_1.fastq',
        test_file_dir / 'assembly' / 'reads_illumina_2.fastq']
    fastq_se_iontorrent = test_file_dir / 'assembly' / 'reads_iontorrent.fastq'
    fastq_se_ont = test_file_dir / 'assembly' / 'ont_bsubtilis_small_region.fastq.gz'

    def test_assembly_illumina(self) -> None:
        """
        Tests the assembly workflow.
        :return: None
        """
        wrapper = AssemblyWrapper(self.running_dir, 'illumina')
        fastq_input = FastqInput('illumina', pe=[ToolIOFile(x) for x in TestWorkflowAssembly.fastq_pe])
        wrapper.run('test_sample', fastq_input, min_ctg_len=None, assembler_opts={'kmers': '25,33'})
        self.assertGreater(wrapper.output.fasta_contigs.stat().st_size, 0)

    def test_assembly_illumina_stats(self) -> None:
        """
        Tests the assembly workflow with standard forward / reverse reads input and with stats determined.
        :return: None
        """
        wrapper = AssemblyWrapper(self.running_dir, 'illumina')
        fastq_input = FastqInput('illumina', pe=[ToolIOFile(x) for x in TestWorkflowAssembly.fastq_pe])
        wrapper.run('test_sample', fastq_input, min_ctg_len=None, assembler_opts={'kmers': '25,33'}, calc_qc_stats=True)
        self.assertGreater(wrapper.output.fasta_contigs.stat().st_size, 0)
        self.assertIsNotNone(wrapper.output.qc_stats)
        self.assertIn('depth', wrapper.output.qc_stats)
        self.assertIn('mapping', wrapper.output.qc_stats)

    def test_assembly_ont(self) -> None:
        """
        Tests the assembly workflow with ONT data.
        :return: None
        """
        wrapper = AssemblyWrapper(self.running_dir, 'ont')
        fastq_input = FastqInput('nanopore', se=[ToolIOFile(TestWorkflowAssembly.fastq_se_ont)], is_pe=False)
        wrapper.run('test_sample', fastq_input, min_ctg_len=500, assembler_opts={'genome_size': '15k'})
        self.assertGreater(wrapper.output.fasta_contigs.stat().st_size, 0)

    def test_assembly_ont_stats(self) -> None:
        """
        Tests the assembly workflow with ONT data and stats generation.
        :return: None
        """
        wrapper = AssemblyWrapper(self.running_dir, 'ont')
        fastq_input = FastqInput('nanopore', se=[ToolIOFile(TestWorkflowAssembly.fastq_se_ont)], is_pe=False)
        wrapper.run('test_sample', fastq_input, min_ctg_len=500, calc_qc_stats=True)
        self.assertGreater(wrapper.output.fasta_contigs.stat().st_size, 0)
        self.assertIsNotNone(wrapper.output.qc_stats)
        self.assertIn('depth', wrapper.output.qc_stats)
        self.assertIn('mapping', wrapper.output.qc_stats)


if __name__ == '__main__':
    unittest.main()
