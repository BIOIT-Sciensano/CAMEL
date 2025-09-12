from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.scripts.mega.mainmega import MainMega


class TestMEGA(CamelTestSuite):
    """
    Tests the MEGA tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('mega')
    input_snp_matrix = test_file_dir / 'sample_snp_matrix.fasta'
    input_vcf_files = [
        test_file_dir / 'variants-s1.vcf',
        test_file_dir / 'variants-s2.vcf',
        test_file_dir / 'variants-s3.vcf',
        test_file_dir / 'variants-s4.vcf'
    ]

    def test_model_selection(self) -> None:
        """
        Tests the model selection
        :return: None
        """
        output_file = self.running_dir / 'out' / 'output_model.tsv'
        if not output_file.parent.exists():
            output_file.parent.mkdir(parents=True)
        args = [
            '--fasta', str(TestMEGA.input_snp_matrix),
            '--working-dir', str(self.running_dir),
            '--action', 'model',
            '--output-model', str(output_file),
            '--missing-data', 'complete_deletion',
            '--site-cov-cutoff', '50',
            '--branch-swap', 'weak'
        ]
        mega = MainMega(args)
        mega.run()
        self.assertGreater(output_file.stat().st_size, 0)

    def test_tree_building(self) -> None:
        """
        Tests the tree building.
        :return: None
        """
        output_file = self.running_dir / 'out' / 'output_tree.nwk'
        if not output_file.parent.exists():
            output_file.parent.mkdir(parents=True)
        args = [
            '--fasta', str(TestMEGA.input_snp_matrix),
            '--working-dir', str(self.running_dir),
            '--action', 'tree',
            '--model', 'T92',
            '--rates', 'G',
            '--output-tree', str(output_file),
            '--missing-data', 'use_all_sites',
            '--branch-swap', 'moderate',
            '--bootstraps', '10',
        ]
        mega = MainMega(args)
        mega.run()
        self.assertGreater(output_file.stat().st_size, 0)

    def test_both(self) -> None:
        """
        Tests model selection + tree building.
        :return: None
        """
        output_file_tree = self.running_dir / 'out' / 'output_tree.nwk'
        output_file_model = self.running_dir / 'out' / 'output_model.tsv'
        if not output_file_model.parent.exists():
            output_file_model.parent.mkdir(parents=True)
        args = [
            '--fasta', str(TestMEGA.input_snp_matrix),
            '--working-dir', str(self.running_dir),
            '--action', 'both',
            '--output-tree', str(output_file_tree),
            '--output-model', str(output_file_model),
            '--bootstraps', '10'
        ]
        mega = MainMega(args)
        mega.run()
        self.assertGreater(output_file_model.stat().st_size, 0)
        self.assertGreater(output_file_tree.stat().st_size, 0)

    def test_export_snp_matrix(self) -> None:
        """
        Tests the export SNP matrix function (starting from multiple VCF input files).
        :return: None
        """
        output_file_model = self.running_dir / 'out' / 'output_model.tsv'
        output_file_fasta = self.running_dir / 'out' / 'snp_matrix.fasta'
        if not output_file_model.parent.exists():
            output_file_model.parent.mkdir(parents=True)
        args = [
            '--output-snp-matrix', str(output_file_fasta),
            '--output-model', str(output_file_model),
            '--working-dir', str(self.running_dir),
            '--action', 'model',
        ]
        for vcf_file in TestMEGA.input_vcf_files:
            args.extend(['--vcf', str(vcf_file), vcf_file.name])
        mega = MainMega(args)
        mega.run()
        self.assertGreater(output_file_model.stat().st_size, 0)
        self.assertGreater(output_file_fasta.stat().st_size, 0)
