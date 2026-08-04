import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import fastqutils

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.bbtools.demuxbyname import Demuxbyname


class TestDemuxbyname(CamelTestSuite):
    """
    Tests the Demuxbyname tool.
    """

    SEGMENTS = ('HA', 'NA', 'NP')
    READS_PER_SEGMENT = 10

    @classmethod
    def _create_test_fastq(cls, path: Path) -> None:
        """
        Creates a synthetic FASTQ file with reads labelled by segment.
        Read names have the format @read_{segment}_{i}/{segment}, so that
        demuxbyname with delimiter='/' and prefixmode=f splits on the segment suffix.
        :param path: Output path for the FASTQ file
        :return: None
        """
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, 'w') as handle:
            for segment in cls.SEGMENTS:
                for i in range(cls.READS_PER_SEGMENT):
                    handle.write(f'@read_{segment}_{i}/{segment}\n')
                    handle.write('ACGTACGTACGTACGTACGT\n')
                    handle.write('+\n')
                    handle.write('IIIIIIIIIIIIIIIIIIII\n')

    def test_demuxbyname(self) -> None:
        """
        Tests Demuxbyname with default parameters (delimiter='/', prefixmode=f).
        Verifies that one output file is created per segment, each containing the
        expected number of reads.
        :return: None
        """
        path_fastq = self.running_dir / 'input' / 'reads_combined.fastq'
        self._create_test_fastq(path_fastq)
        demuxbyname = Demuxbyname()
        demuxbyname.add_input_files({'FASTQ': [ToolIOFile(path_fastq)]})
        demuxbyname.run(self.running_dir)
        print(demuxbyname.tool_outputs)
        self.assertEqual(len(demuxbyname.tool_outputs['FASTQ']), len(TestDemuxbyname.SEGMENTS))
        for fastq_out in demuxbyname.tool_outputs['FASTQ']:
            self.assertEqual(fastqutils.count_reads(fastq_out.path), TestDemuxbyname.READS_PER_SEGMENT)

    def test_demuxbyname_read_count(self) -> None:
        """
        Tests that the total number of reads across all output files equals the
        number of reads in the input file.
        :return: None
        """
        path_fastq = self.running_dir / 'input' / 'reads_combined.fastq'
        self._create_test_fastq(path_fastq)
        nb_reads_in = fastqutils.count_reads(path_fastq)
        demuxbyname = Demuxbyname()
        demuxbyname.add_input_files({'FASTQ': [ToolIOFile(path_fastq)]})
        demuxbyname.run(self.running_dir)
        nb_reads_out = sum(fastqutils.count_reads(f.path) for f in demuxbyname.tool_outputs['FASTQ'])
        self.assertEqual(nb_reads_out, nb_reads_in)

    def test_demuxbyname_no_input(self) -> None:
        """
        Tests that Demuxbyname raises an error when no FASTQ input is provided.
        :return: None
        """
        demuxbyname = Demuxbyname()
        with self.assertRaises(InvalidToolInputError):
            demuxbyname.run(self.running_dir)


if __name__ == '__main__':
    unittest.main()
