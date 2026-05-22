import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurdistseqs import MothurDistSeqs


class TestDistSeqs(CamelTestSuite):
    """
    Tests Mothur dist.seqs.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_FASTA = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta')

    def test_distseqs(self) -> None:
        """
        Tests Mothur dist.seqs.
        :return: None
        """
        distseqs = MothurDistSeqs()
        distseqs.add_input_files({
            'FASTA': [TestDistSeqs.FILE_FASTA]
        })
        distseqs.run(self.running_dir)
        self.assertTrue('DIST' in distseqs.tool_outputs, "No distance matrix output generated")
        chimera_output = Path(distseqs.tool_outputs['DIST'][0].path)
        self.assertTrue(chimera_output.exists())
        self.assertGreater(chimera_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
