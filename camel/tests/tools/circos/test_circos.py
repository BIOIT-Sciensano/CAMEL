import unittest

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.circos.circos import Circos


class TestCircos(CamelTestSuite):
    """
    Tests Circos.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('circos')
    FILE_CIRCOS_CONFIG = ToolIOFile(test_file_dir / 'hello_world.txt')

    def test_circos(self) -> None:
        """
        Tests circos.
        :return: None
        """
        circos = Circos()
        circos.add_input_files({'TXT': [TestCircos.FILE_CIRCOS_CONFIG]})
        circos.run(self.running_dir)
        self.assertTrue(len(circos.tool_outputs) > 0, "No outputs generated")


if __name__ == '__main__':
    unittest.main()
