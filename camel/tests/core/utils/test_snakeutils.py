import unittest
from importlib.resources import files
from pathlib import Path

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.snakemake import snakepipelineutils


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


if __name__ == '__main__':
    unittest.main()
