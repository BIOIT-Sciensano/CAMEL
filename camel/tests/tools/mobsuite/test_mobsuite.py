import unittest
from pathlib import Path

from camelcore.app.io.tooliodirectory import ToolIODirectory
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.config import config
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mobsuite.mobrecon import MOBRecon
from camel.app.tools.mobsuite.mobreconreporter import MOBReconReporter


class TestMOBSuite(CamelTestSuite):
    """
    Tests the MOB-suite tool.
    """
    # Get test file and reference file directories
    test_file_dir = CamelTestSuite.get_test_file_dir('mob_suite')
    input_fasta_plasmid = test_file_dir / 'AB011548.fasta'
    input_fasta_contigs = test_file_dir / 'ecoli_contigs.fasta'
    input_fasta_no_plasmid = test_file_dir / 'SRR15000905-ds_contigs.fasta'
    input_fasta_multi_contigs = test_file_dir / 'multi_contigs.fasta'

    def test_mob_recon(self) -> None:
        """
        Tests the MOB-recon tool.
        :return: None
        """
        mob_recon = MOBRecon()
        mob_recon.add_input_files({
            'FASTA': [ToolIOFile(TestMOBSuite.input_fasta_plasmid)],
            'DB': [ToolIODirectory(Path(config.dir_db / 'mob_suite' / 'latest'))]
        })
        mob_recon.update_parameters(num_threads=8)
        mob_recon.run(self.running_dir)
        self.verify_output_files(mob_recon, 'TSV')
        self.verify_output_files(mob_recon, 'FASTA')
        self.assertIn('detected_plasmids', mob_recon.informs)

    def test_mob_recon_assembly(self) -> None:
        """
        Tests the MOB-recon tool on assembled contigs.
        :return: None
        """
        mob_recon = MOBRecon()
        mob_recon.add_input_files({
            'FASTA': [ToolIOFile(TestMOBSuite.input_fasta_contigs)],
            'DB': [ToolIODirectory(Path(config.dir_db / 'mob_suite' / 'latest'))]
        })
        mob_recon.update_parameters(num_threads=8)
        mob_recon.run(self.running_dir)
        self.verify_output_files(mob_recon, 'TSV')
        self.verify_output_files(mob_recon, 'FASTA', nb_files=2)
        self.assertIn('detected_plasmids', mob_recon.informs)

    def test_mob_recon_reporter(self) -> None:
        """
        Tests the reporter for the MOB-recon tool.
        :return: None
        """
        # Run MOB-recon
        mob_recon = MOBRecon()
        mob_recon.add_input_files({
            'FASTA': [ToolIOFile(TestMOBSuite.input_fasta_plasmid)],
            'DB': [ToolIODirectory(Path(config.dir_db / 'mob_suite' / 'latest'))]
        })
        mob_recon.update_parameters(num_threads=8)
        mob_recon.run(self.running_dir)

        # Run the reporter
        reporter = MOBReconReporter()
        reporter.add_input_files({
            'TSV': mob_recon.tool_outputs['TSV'],
            'TSV_contigs': mob_recon.tool_outputs['TSV_contigs'],
            'FASTA': mob_recon.tool_outputs['FASTA']})
        reporter.add_input_informs({'mob_recon': mob_recon.informs})
        reporter.run(self.running_dir)

        # Check the output
        output_section = reporter.tool_outputs['HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)
        self.export_report_section(reporter.tool_outputs['HTML'][0].value, self.running_dir / 'report')

    def test_mob_recon_reporter_with_contig_report(self) -> None:
        """
        Tests the reporter for the MOB-recon tool.
        :return: None
        """
        # Run MOB-recon
        mob_recon = MOBRecon()
        mob_recon.add_input_files({
            'FASTA': [ToolIOFile(TestMOBSuite.input_fasta_multi_contigs)],
            'DB': [ToolIODirectory(Path(config.dir_db / 'mob_suite' / 'latest'))]
        })
        mob_recon.update_parameters(num_threads=8)
        mob_recon.run(self.running_dir)

        # Run the reporter
        reporter = MOBReconReporter()
        reporter.add_input_files({
            'TSV': mob_recon.tool_outputs['TSV'],
            'TSV_contigs': mob_recon.tool_outputs['TSV_contigs'],
            'FASTA': mob_recon.tool_outputs['FASTA']})
        reporter.add_input_informs({'mob_recon': mob_recon.informs})
        reporter.update_parameters(contig_report=True)
        reporter.run(self.running_dir)

        # Check the output
        output_section = reporter.tool_outputs['HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)
        self.export_report_section(reporter.tool_outputs['HTML'][0].value, self.running_dir / 'report')

    def test_mob_recon_assembly_no_plasmid(self) -> None:
        """
        Tests the MOB-recon tool on assembled contigs without any plasmids.
        :return: None
        """
        mob_recon = MOBRecon()
        mob_recon.add_input_files({
            'FASTA': [ToolIOFile(TestMOBSuite.input_fasta_no_plasmid)],
            'DB': [ToolIODirectory(Path(config.dir_db / 'mob_suite' / 'latest'))]
        })
        mob_recon.update_parameters(num_threads=8)
        mob_recon.run(self.running_dir)
        # Empty file -> verify_output_files cannot be used
        self.assertTrue(mob_recon.tool_outputs['TSV'][0].path.exists())
        self.verify_output_files(mob_recon, 'FASTA', nb_files=0)

    def test_mob_recon_reporter_no_plasmid(self) -> None:
        """
        Tests the reporter for the MOB-recon tool on contigs without any plasmids.
        :return: None
        """
        # Run MOB-recon
        mob_recon = MOBRecon()
        mob_recon.add_input_files({
            'FASTA': [ToolIOFile(TestMOBSuite.input_fasta_no_plasmid)],
            'DB': [ToolIODirectory(Path(config.dir_db / 'mob_suite' / 'latest'))]
        })
        mob_recon.update_parameters(num_threads=8)
        mob_recon.run(self.running_dir)

        # Run the reporter
        reporter = MOBReconReporter()
        reporter.add_input_files({
            'TSV': mob_recon.tool_outputs['TSV'],
            'TSV_contigs': mob_recon.tool_outputs['TSV_contigs'],
            'FASTA': mob_recon.tool_outputs['FASTA']})
        reporter.add_input_informs({'mob_recon': mob_recon.informs})
        reporter.run(self.running_dir)

        # Check the output
        output_section = reporter.tool_outputs['HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
