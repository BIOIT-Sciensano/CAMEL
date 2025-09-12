from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.components.workflows.genedetectionwrapper import GeneDetectionWrapper
from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.io.tooliofile import ToolIOFile


class TestWorkflowGeneDetection(CamelTestSuite):
    """
    Tests the Illumina trimming workflow.
    """
    # Input files (Gene detection)
    test_file_dir = CamelTestSuite.get_test_file_dir('gene_detection')

    # Input files
    gene_db = test_file_dir / 'db'
    input_fastq_by_key = {
        'illumina': [test_file_dir / 'illumina' / 'reads_illumina_1.fastq',
                     test_file_dir / 'illumina' / 'reads_illumina_2.fastq'],
        'iontorrent': [test_file_dir / 'iontorrent' / 'reads_iontorrent.fastq'],
        'nanopore': [test_file_dir / 'nanopore' / 'reads_nanopore.fastq']
    }
    fasta = test_file_dir / 'contigs.fasta'

    def test_gene_detection_blast(self) -> None:
        """
        Tests the gene detection workflow with BLAST+ detection and FASTA input.
        :return: None
        """
        wrapper = GeneDetectionWrapper(self.running_dir)
        db_data = {
            'path': str(TestWorkflowGeneDetection.gene_db),
            'min_percent_identity': 80,
            'min_coverage': 80
        }
        wrapper.run_workflow_blast(TestWorkflowGeneDetection.fasta, 'test sample', db_data)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)

    def test_gene_detection_illumina_kma(self) -> None:
        """
        Tests the gene detection workflow with KMA detection and Illumina data.
        :return: None
        """
        wrapper = GeneDetectionWrapper(self.running_dir)
        fastq_input = FastqInput('illumina', pe=[
            ToolIOFile(x) for x in TestWorkflowGeneDetection.input_fastq_by_key['illumina']])
        db_data = {
            'path': str(TestWorkflowGeneDetection.gene_db),
            'min_percent_identity': 80,
            'min_coverage': 80
        }
        wrapper.run_workflow_kma(fastq_input, 'test sample', db_data)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)

    def test_gene_detection_nanopore_kma(self) -> None:
        """
        Tests the gene detection workflow with KMA detection and IonTorrent data.
        :return: None
        """
        wrapper = GeneDetectionWrapper(self.running_dir)
        fastq_input = FastqInput('ont', se=[
            ToolIOFile(TestWorkflowGeneDetection.input_fastq_by_key['nanopore'][0])], is_pe=False)
        db_data = {
            'path': str(TestWorkflowGeneDetection.gene_db),
            'min_percent_identity': 80,
            'min_coverage': 80
        }
        wrapper.run_workflow_kma(fastq_input, 'test sample', db_data)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)
