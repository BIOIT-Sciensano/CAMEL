import filecmp
import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.errors import ToolExecutionError, InvalidParameterError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.muscle.muscle import Muscle


class TestMuscle(CamelTestSuite):
    """
    Tests the Muscle tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('muscle')

    def test_muscle(self) -> None:
        """
        Tests basic Muscle run.
        :return: None
        """
        muscle = Muscle()
        muscle.add_input_files({
            'FASTA': [ToolIOFile(TestMuscle.test_file_dir / "unaligned.fasta")]
        })
        muscle.run(self.running_dir)
        self.verify_output_files(muscle, 'FASTA')

        compfile = TestMuscle.test_file_dir / 'aligned.fasta'
        self.assertTrue(filecmp.cmp(muscle.tool_outputs['FASTA'][0].path, compfile, shallow=False))

    def test_muscle_html(self) -> None:
        """
        Tests basic Muscle run with html output.
        :return: None
        """
        muscle = Muscle()
        muscle.add_input_files({
            'FASTA': [ToolIOFile(TestMuscle.test_file_dir / "unaligned.fasta")]
        })
        muscle.update_parameters(html=None)
        muscle.run(self.running_dir)
        self.verify_output_files(muscle, 'HTML')

        compfile = TestMuscle.test_file_dir / 'aligned.html'
        self.assertTrue(filecmp.cmp(muscle.tool_outputs['HTML'][0].path, compfile, shallow=False))

    def test_muscle_clw(self) -> None:
        """
        Tests basic Muscle run with clw output.
        :return: None
        """
        muscle = Muscle()
        muscle.add_input_files({
            'FASTA': [ToolIOFile(TestMuscle.test_file_dir / "unaligned.fasta")]
        })
        muscle.update_parameters(clw=None)
        muscle.run(self.running_dir)
        self.verify_output_files(muscle, 'CLW')

        compfile = TestMuscle.test_file_dir / 'aligned.clw'
        self.assertTrue(filecmp.cmp(muscle.tool_outputs['CLW'][0].path, compfile, shallow=False))

    def test_muscle_clw_strict(self) -> None:
        """
        Tests basic Muscle run with clw strict output.
        :return: None
        """
        muscle = Muscle()
        muscle.add_input_files({
            'FASTA': [ToolIOFile(TestMuscle.test_file_dir / "unaligned.fasta")]
        })
        muscle.update_parameters(clwstrict=None)
        muscle.run(self.running_dir)
        self.verify_output_files(muscle, 'CLW')

        compfile = TestMuscle.test_file_dir / 'aligned.clws'
        self.assertTrue(filecmp.cmp(muscle.tool_outputs['CLW'][0].path, compfile, shallow=False))

    def test_muscle_msf(self) -> None:
        """
        Tests basic Muscle run with msf output.
        :return: None
        """
        muscle = Muscle()
        muscle.add_input_files({
            'FASTA': [ToolIOFile(TestMuscle.test_file_dir / "unaligned.fasta")]
        })
        muscle.update_parameters(msf=None)
        muscle.run(self.running_dir)
        self.verify_output_files(muscle, 'MSF')

        compfile = TestMuscle.test_file_dir / 'aligned.msf'
        self.assertTrue(filecmp.cmp(muscle.tool_outputs['MSF'][0].path, compfile, shallow=False))

    def test_muscle_error(self) -> None:
        """
        Tests Muscle with broken fasta file.
        :return: None
        """
        muscle = Muscle()
        muscle.add_input_files({
            'FASTA': [ToolIOFile(TestMuscle.test_file_dir / "bad.fasta")]
        })

        with self.assertRaises(ToolExecutionError):
            muscle.run(self.running_dir)

    def test_muscle_log_output(self) -> None:
        """
        Tests Muscle run with log output.
        :return: None
        """
        muscle = Muscle()
        muscle.add_input_files({
            'FASTA': [ToolIOFile(TestMuscle.test_file_dir / "unaligned.fasta")]
        })
        muscle.update_parameters(log=self.running_dir / 'log.txt')
        muscle.run(self.running_dir)
        self.verify_output_files(muscle, 'LOG')

    def test_muscle_loga_output(self) -> None:
        """
        Tests Muscle run with appending the log output.
        :return: None
        """
        muscle = Muscle()
        muscle.add_input_files({
            'FASTA': [ToolIOFile(TestMuscle.test_file_dir / "unaligned.fasta")]
        })
        muscle.update_parameters(loga=self.running_dir / 'log.txt')
        muscle.run(self.running_dir)
        self.verify_output_files(muscle, 'LOG')

    def test_muscle_output_param_error(self) -> None:
        """
        Tests Muscle with too many output formats specified.
        :return: None
        """
        muscle = Muscle()
        muscle.add_input_files({
            'FASTA': [ToolIOFile(TestMuscle.test_file_dir / "unaligned.fasta")]
        })
        muscle.update_parameters(msf=None, html=None)

        with self.assertRaises(InvalidParameterError):
            muscle.run(self.running_dir)

    def test_muscle_log_param_error(self) -> None:
        """
        Tests Muscle with too many log flags specified.
        :return: None
        """
        muscle = Muscle()
        muscle.add_input_files({
            'FASTA': [ToolIOFile(TestMuscle.test_file_dir / "unaligned.fasta")]
        })
        muscle.update_parameters(log=self.running_dir / 'log.txt', loga=self.running_dir / 'loga.txt')

        with self.assertRaises(InvalidParameterError):
            muscle.run(self.running_dir)


if __name__ == '__main__':
    unittest.main()
