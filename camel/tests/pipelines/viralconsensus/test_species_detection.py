import unittest
from pathlib import Path

from Bio import SeqIO
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.config import config
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.scriptutils.basepipe.fastqinput import FastqInput
from camel.scripts.viralconsensuspipeline.workflows import speciesdetection
from camel.scripts.viralconsensuspipeline.workflows.speciesdetection import SpeciesDetectionWorkflow
from camel.tests import longRunningTest


class TestSpeciesDetection(CamelTestSuite):
    """
    Tests the species detection for the viral consensus pipeline.
    """

    dir_testdata = CamelTestSuite.get_test_file_dir('pipelines', 'viral_consensus')
    dir_ref_genomes = Path(config.dir_db, 'pipelines', 'viral_consensus', 'version_1.1', 'ref_genomes')
    SPECIES_BY_SEQ_ID = {'refA1-genome': 'rsv_a', 'refA2-genome': 'rsv_a', 'refB-genome': 'rsv_b'}

    def _create_depth_file(self, depth_by_seq_id: dict[str, list[int]]) -> Path:
        """
        Creates a samtools depth file.
        :param depth_by_seq_id: Depth values by sequence id
        :return: Path to the depth file
        """
        path_depth = self.running_dir / 'depth.tsv'
        with path_depth.open('w') as handle:
            for seq_id, depths in depth_by_seq_id.items():
                for pos, depth in enumerate(depths, start=1):
                    handle.write(f'{seq_id}\t{pos}\t{depth}\n')
        return path_depth

    def test_create_reference_fasta(self) -> None:
        """
        Tests combining the reference genomes of the candidate species.
        :return: None
        """
        path_out = self.running_dir / 'references.fasta'
        species_by_seq_id = speciesdetection.create_reference_fasta({
            'rsv_a': TestSpeciesDetection.dir_ref_genomes / 'hRSV_A.fasta',
            'rsv_b': TestSpeciesDetection.dir_ref_genomes / 'hRSV_B.fasta',
        }, path_out)
        self.assertEqual(sorted(species_by_seq_id.values()), ['rsv_a', 'rsv_b'])
        with path_out.open() as handle:
            self.assertEqual(len(list(SeqIO.parse(handle, 'fasta'))), 2)

    def test_summarize_depth(self) -> None:
        """
        Tests summarizing the depth by species.
        :return: None
        """
        path_depth = self._create_depth_file({
            'refA1-genome': [10, 10, 0, 0],
            'refA2-genome': [5, 5, 5, 5],
            'refB-genome': [0, 0, 0, 0],
        })
        data_species = speciesdetection.summarize_depth(path_depth, TestSpeciesDetection.SPECIES_BY_SEQ_ID)
        self.assertEqual(list(data_species['species']), ['rsv_a', 'rsv_b'])
        self.assertEqual(list(data_species['mapped_bases']), [40, 0])
        self.assertEqual(list(data_species['fraction_bases']), [1.0, 0.0])
        self.assertEqual(list(data_species['covered_rate']), [1.0, 0.0])

    def test_select_species(self) -> None:
        """
        Tests selecting the species with the most mapped bases.
        :return: None
        """
        path_depth = self._create_depth_file({
            'refA1-genome': [1, 1, 0, 0],
            'refA2-genome': [0, 0, 0, 0],
            'refB-genome': [50, 50, 50, 50],
        })
        data_species = speciesdetection.summarize_depth(path_depth, TestSpeciesDetection.SPECIES_BY_SEQ_ID)
        self.assertEqual(speciesdetection.select_species(data_species), 'rsv_b')

    def test_select_species_ambiguous(self) -> None:
        """
        Tests that no species is selected when the mapped bases are split between species.
        :return: None
        """
        path_depth = self._create_depth_file({
            'refA1-genome': [10, 10, 10, 10],
            'refA2-genome': [0, 0, 0, 0],
            'refB-genome': [8, 8, 8, 8],
        })
        data_species = speciesdetection.summarize_depth(path_depth, TestSpeciesDetection.SPECIES_BY_SEQ_ID)
        self.assertIsNone(speciesdetection.select_species(data_species))

    def test_select_species_no_reads(self) -> None:
        """
        Tests that no species is selected when no reads are mapped.
        :return: None
        """
        path_depth = self._create_depth_file({
            'refA1-genome': [0, 0, 0, 0],
            'refA2-genome': [0, 0, 0, 0],
            'refB-genome': [0, 0, 0, 0],
        })
        data_species = speciesdetection.summarize_depth(path_depth, TestSpeciesDetection.SPECIES_BY_SEQ_ID)
        self.assertIsNone(speciesdetection.select_species(data_species))

    @longRunningTest()
    def test_species_detection_workflow_illumina(self) -> None:
        """
        Tests the species detection workflow with Illumina RSV A data.
        :return: None
        """
        fastq_in = FastqInput('illumina', pe=[
            ToolIOFile(TestSpeciesDetection.dir_testdata / 'RSV-ERR331022_1.fastq.gz'),
            ToolIOFile(TestSpeciesDetection.dir_testdata / 'RSV-ERR331022_2.fastq.gz'),
        ], is_pe=True)
        workflow = SpeciesDetectionWorkflow(self.running_dir)
        out = workflow.run(fastq_in, {
            'rsv_a': TestSpeciesDetection.dir_ref_genomes / 'hRSV_A.fasta',
            'rsv_b': TestSpeciesDetection.dir_ref_genomes / 'hRSV_B.fasta',
        })
        self.assertEqual(out.species, 'rsv_a')
        self.assertTrue(out.path_tsv.exists())
        self.assertTrue(out.path_informs.exists())


if __name__ == '__main__':
    unittest.main()
