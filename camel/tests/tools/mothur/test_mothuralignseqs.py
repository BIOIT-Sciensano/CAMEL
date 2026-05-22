import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothuralignseqs import MothurAlignSeqs


class TestMothurAlignSeqs(CamelTestSuite):
    """
    Tests the MothurAlignSeqs tool.
    """

    def test_init(self) -> None:
        """
        Tests that the tool initializes correctly and the version is retrieved.
        This is a placeholder until real tests are written.
        """
        tool = MothurAlignSeqs()
        self.assertIsNotNone(tool.version)


if __name__ == '__main__':
    unittest.main()
