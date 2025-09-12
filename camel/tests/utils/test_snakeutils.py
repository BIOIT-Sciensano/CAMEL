import unittest
from importlib.resources import files
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils


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
        config_path = SnakePipelineUtils.generate_config_file(config_data, self.running_dir)
        SnakePipelineUtils.run_snakemake(str(TestSnakeUtils.test_snakefile), config_path, [], self.running_dir)

    def test_run_snakefile_resources(self) -> None:
        """
        Tests the run snakefile method with the resources parameter set.
        :return: None
        """
        config_data = {
            'working_dir': str(self.running_dir)
        }
        config_path = SnakePipelineUtils.generate_config_file(config_data, self.running_dir)
        command = SnakePipelineUtils.run_snakemake(
            str(TestSnakeUtils.test_snakefile), config_path, [], self.running_dir, resources={'RAM': 2, 'GPU': 4})
        self.assertIn('--resources', command.command, 'Resources parameter not added')
        self.assertIn('RAM', command.command, 'Resource not added')
        self.assertIn('GPU', command.command, 'Resource not added')


if __name__ == '__main__':
    unittest.main()
