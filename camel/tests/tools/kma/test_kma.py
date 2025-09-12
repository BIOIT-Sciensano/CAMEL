import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.kma.kma import KMA


class TestKMA(CamelTestSuite):
    """
    Tests the KMA tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('kma')
    input_reads = {
        'illumina': [test_file_dir / 'reads_illumina_1.fastq', test_file_dir / 'reads_illumina_2.fastq'],
        'iontorrent': [test_file_dir / 'reads_iontorrent.fastq'],
        'nanopore': [test_file_dir / 'reads_nanopore.fastq']
    }
    input_gene_detection_db = test_file_dir / 'db' / 'arg-annot-clustered_80'

    def test_kma_pe(self) -> None:
        """
        Tests KMA with paired-end input.
        :return: None
        """
        kma = KMA()
        kma.add_input_files({
            'FASTQ_PE': [ToolIOFile(file_) for file_ in TestKMA.input_reads['illumina']],
            'DB': [ToolIOValue(str(TestKMA.input_gene_detection_db))]
        })
        kma.run(self.running_dir)
        self.assertTrue('TSV' in kma.tool_outputs, "No TSV output generated")
        output_file = Path(kma.tool_outputs['TSV'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_kma_se(self) -> None:
        """
        Tests KMA with single-end input.
        :return: None
        """
        kma = KMA()
        kma.add_input_files({
            'FASTQ_SE': [ToolIOFile(file_) for file_ in TestKMA.input_reads['iontorrent']],
            'DB': [ToolIOValue(str(TestKMA.input_gene_detection_db))]
        })
        kma.run(self.running_dir)
        self.assertTrue('TSV' in kma.tool_outputs, "No TSV output generated")
        output_file = Path(kma.tool_outputs['TSV'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_kma_se_nanopore(self) -> None:
        """
        Tests KMA with single-end nanopore input.
        :return: None
        """
        kma = KMA()
        kma.add_input_files({
            'FASTQ_SE': [ToolIOFile(file_) for file_ in TestKMA.input_reads['nanopore']],
            'DB': [ToolIOValue(str(TestKMA.input_gene_detection_db))]
        })
        kma.update_parameters(bc_nano=None, basecalls='0.7')
        kma.run(self.running_dir)
        self.assertTrue('TSV' in kma.tool_outputs, "No TSV output generated")
        output_file = Path(kma.tool_outputs['TSV'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_kma_ont_preset(self) -> None:
        """
        Tests KMA with preset.
        :return: None
        """
        kma = KMA()
        kma.add_input_files({
            'FASTQ_SE': [ToolIOFile(file_) for file_ in TestKMA.input_reads['nanopore']],
            'DB': [ToolIOValue(str(TestKMA.input_gene_detection_db))]
        })
        kma.update_parameters(bc_nano=None, basecalls='0.7', ont=None)
        kma.run(self.running_dir)
        self.verify_output_files(kma, 'TSV')

    def test_kma_apm_preset(self) -> None:
        """
        Tests KMA with paired-end input.
        :return: None
        """
        kma = KMA()
        kma.add_input_files({
            'FASTQ_PE': [ToolIOFile(file_) for file_ in TestKMA.input_reads['illumina']],
            'DB': [ToolIOValue(str(TestKMA.input_gene_detection_db))]
        })
        kma.update_parameters(apm='p')
        kma.run(self.running_dir)
        self.verify_output_files(kma, 'TSV')

if __name__ == '__main__':
    unittest.main()
