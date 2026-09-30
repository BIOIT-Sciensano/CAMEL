import dataclasses
import json
import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from camel.app.scriptutils.basepipe.fastqinput import FastqInput
from camel.scripts.viralconsensuspipeline.workflows.readmappingworkflow import ReadMappingWorkflow

# Minimum fraction of the mapped bases that should map to the detected species
MIN_FRACTION_BASES = 0.8

# Minimum fraction of the reference genome(s) of the detected species that should be covered by at least one read
MIN_COVERED_RATE = 0.05


@dataclasses.dataclass
class SpeciesDetectionOutput:
    """
    Holder for the output of the species detection workflow.
    """
    species: str | None
    path_tsv: Path
    path_informs: Path


class SpeciesDetectionWorkflow:
    """
    Detects the species by competitive read mapping against a reference genome for each candidate species.
    """

    def __init__(self, dir_: Path) -> None:
        """
        Initializes this workflow.
        :param dir_: Working directory
        :return: None
        """
        self._dir = dir_
        if not self._dir.exists():
            logging.info(f'Creating working directory: {self._dir}')
            self._dir.mkdir(parents=True)

    def run(self, fastq_in: FastqInput, fasta_by_species: dict[str, Path], threads: int = 8) -> SpeciesDetectionOutput:
        """
        Runs the species detection workflow.
        :param fastq_in: Input FASTQ data
        :param fasta_by_species: Reference genome FASTA file by candidate species
        :param threads: Number of threads
        :return: Output holder
        """
        # Map the reads against the combined references
        path_fasta = self._dir / 'references.fasta'
        species_by_seq_id = create_reference_fasta(fasta_by_species, path_fasta)
        out_mapping = ReadMappingWorkflow(self._dir / 'mapping').run(fastq_in, path_fasta, threads=threads)

        # Select the species
        data_species = summarize_depth(out_mapping.path_depth, species_by_seq_id)
        species = select_species(data_species)

        # Save output
        path_tsv = self._dir / 'species_detection.tsv'
        data_species.to_csv(path_tsv, sep='\t', index=False)
        path_informs = self._dir / 'informs.json'
        with path_informs.open('w') as handle:
            json.dump(out_mapping.informs, handle, indent=2)
        return SpeciesDetectionOutput(species=species, path_tsv=path_tsv, path_informs=path_informs)


def create_reference_fasta(fasta_by_species: dict[str, Path], path_out: Path) -> dict[str, str]:
    """
    Combines the reference genomes of the candidate species into a single FASTA file.
    :param fasta_by_species: Reference genome FASTA file by candidate species
    :param path_out: Output FASTA path
    :return: Species by sequence id
    """
    species_by_seq_id = {}
    records = []
    for species, path_fasta in fasta_by_species.items():
        with path_fasta.open() as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                if record.id in species_by_seq_id:
                    raise ValueError(f'Duplicate sequence id in the species detection references: {record.id}')
                species_by_seq_id[record.id] = species
                records.append(record)
    with path_out.open('w') as handle:
        SeqIO.write(records, handle, 'fasta')
    return species_by_seq_id


def summarize_depth(path_depth: Path, species_by_seq_id: dict[str, str]) -> pd.DataFrame:
    """
    Summarizes the read depth by species.
    :param path_depth: Samtools depth output (all positions)
    :param species_by_seq_id: Species by sequence id
    :return: Mapped bases, fraction of the mapped bases and covered rate by species (sorted on the mapped bases)
    """
    data_depth = pd.read_table(path_depth, names=['seq_id', 'pos', 'depth'])
    data_by_seq = data_depth.groupby('seq_id')['depth'].agg(
        mapped_bases='sum', covered_rate=lambda depth: (depth > 0).mean())
    data_by_seq['species'] = data_by_seq.index.map(species_by_seq_id)
    data_species = data_by_seq.groupby('species').agg(
        mapped_bases=('mapped_bases', 'sum'), covered_rate=('covered_rate', 'max'))
    data_species = data_species.reindex(sorted(set(species_by_seq_id.values())), fill_value=0)
    total_bases = data_species['mapped_bases'].sum()
    data_species['fraction_bases'] = data_species['mapped_bases'] / total_bases if total_bases > 0 else 0.0
    data_species = data_species.rename_axis('species').reset_index()
    return data_species.sort_values(by='mapped_bases', ascending=False, ignore_index=True)


def select_species(data_species: pd.DataFrame) -> str | None:
    """
    Selects the species with the most mapped bases, if it passes the detection cutoffs.
    :param data_species: Species depth summary (sorted on the mapped bases)
    :return: Detected species, None if no species could be detected
    """
    best = data_species.iloc[0]
    if best['covered_rate'] < MIN_COVERED_RATE:
        logging.warning(f"Insufficient coverage for {best['species']} ({best['covered_rate']:.2%})")
        return None
    if best['fraction_bases'] < MIN_FRACTION_BASES:
        logging.warning(f"Ambiguous species detection: {best['fraction_bases']:.2%} of bases on {best['species']}")
        return None
    return str(best['species'])
