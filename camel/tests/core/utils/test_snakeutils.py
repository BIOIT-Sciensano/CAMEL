import unittest
from importlib.resources import files
from pathlib import Path

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.snakemake import snakepipelineutils, snakemakeutils


class TestSnakeUtils(CamelTestSuite):
    """
    Tests the snakemake utility classes.
    """
    # Input files
    test_snakefile = Path(str(files('camel').joinpath('resources/testing/workflow_test.smk')))

    def test_run_snakefile(self) -> None:
        """
        Tests the run snakefile method.
        :return: None
        """
        config_data = {
            'working_dir': str(self.running_dir)
        }
        config_path = snakepipelineutils.generate_config_file(config_data, self.running_dir)
        snakepipelineutils.run_snakemake(str(TestSnakeUtils.test_snakefile), config_path, [], self.running_dir)

    def test_run_snakefile_resources(self) -> None:
        """
        Tests the run snakefile method with the resources parameter set.
        :return: None
        """
        config_data = {
            'working_dir': str(self.running_dir)
        }
        config_path = snakepipelineutils.generate_config_file(config_data, self.running_dir)
        command = snakepipelineutils.run_snakemake(
            str(TestSnakeUtils.test_snakefile), config_path, [], self.running_dir, resources={'RAM': 2, 'GPU': 4})
        self.assertIn('--resources', command.command, 'Resources parameter not added')
        self.assertIn('RAM', command.command, 'Resource not added')
        self.assertIn('GPU', command.command, 'Resource not added')


    def test_dump_and_load_tooliofile(self) -> None:
        """
        Tests dumping and loading a ToolIOFile object and checks for content/hash preservation.
        :return: None
        """
        test_file_path = CamelTestSuite.get_test_file_dir('components').joinpath('toy.fasta')
        original_io_file = ToolIOFile(path=test_file_path)

        # Save in an IO file
        path_io = self.running_dir / 'fasta.io'
        snakemakeutils.dump_object([original_io_file], path_io)

        # Load from the IO file
        objs_out = snakemakeutils.load_object(path_io)
        loaded_io_file = objs_out[0]

        # Verify that the files match
        self.assertIsInstance(loaded_io_file, ToolIOFile)
        self.assertEqual(loaded_io_file.path.name, original_io_file.path.name)
        self.assertEqual(original_io_file.hash, loaded_io_file.hash)


if __name__ == '__main__':
    unittest.main()
