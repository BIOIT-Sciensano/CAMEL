import argparse
import itertools
import re
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional, Union

import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable, IUPACData
from Bio.Seq import Seq

from camel.app.loggers import logger, initialize_logging
from camel.app.toolkits.mycobacterium.amrutils import ConfidenceLevel


# Example usage:
# update_amr_db.py \
#   --version '2023.6'
#   --tsv-coords /db/pipelines/mycobacterium/amr/who_2023.6/input/mutation_locations_who.tsv
#   --tsv-catalogue /db/pipelines/mycobacterium/amr/who_2023.6/input/amr_associations_who.tsv
#   --tsv-inhouse /db/pipelines/mycobacterium/amr/who_2023.6/input/amr_associations_inhouse.tsv
#   --tsv-abs /db/pipelines/mycobacterium/amr/who_2023.6/input/antibiotics.tsv
#   --ref-fasta /db/refgenomes/Mycobacterium_tuberculosis/H37Rv.fasta
#   --ref-gff /db/refgenomes/Mycobacterium_tuberculosis/H37Rv.gff

class UpdateAMRDB:
    """
    Class to update the Mycobacterium AMR database.
    """

    MAPPING_CONFIDENCE = {
        '1) Assoc w R': ConfidenceLevel.ASSOC_R,
        '2) Assoc w R - Interim': ConfidenceLevel.ASSOC_R_int,
        '3) Uncertain significance': ConfidenceLevel.UNKNOWN,
        '4) Not assoc w R - Interim': ConfidenceLevel.ASSOC_S_int,
        '5) Not assoc w R': ConfidenceLevel.ASSOC_S
    }

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the updater.
        :param args: Command line arguments
        :return: None
        """
        self._args = UpdateAMRDB._parse_arguments(args)
        self._data_mut_locations = None
        self._data_amr_catalogue = None
        self._data_amr_inhouse = None
        self._data_antibiotics = None
        self._seq_ref = None
        self._gff_by_gene_name = {}
        self._gff_by_locus_tag = {}

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments (optional)
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        # Input
        parser.add_argument(
            '--tsv-coords', required=True, type=Path,
            help='TSV file with the genomic coordinates of the mutations')
        parser.add_argument(
            '--tsv-catalogue', required=True, type=Path,
            help='TSV file with the WHO mutation catalogue')
        parser.add_argument(
            '--tsv-inhouse', required=True, type=Path,
            help='TSV file with the in-house mutation catalogue')
        parser.add_argument(
            '--tsv-abs', required=True, type=Path,
            help='TSV file with antibiotics')

        # Reference genome
        parser.add_argument(
            '--ref-gff', required=True, type=Path,
            help='GFF3 file with the ref. genome annotations, used to extract regions')
        parser.add_argument(
            '--ref-fasta', required=True, type=Path, help='Ref. genome FASTA file')
        parser.add_argument('--version', type=str, help='DB version', required=True)

        # Other parameters
        parser.add_argument('--dir-out', type=Path, help='Output directory')
        return parser.parse_args(args)

    @staticmethod
    def get_mutations(codon_in: str, codon_targets: list[str]) -> list:
        """
        Returns all possible mutations that turn the input codon into any of the target codons.
        :param codon_in: Input codon
        :param codon_targets: List of target codons
        :return: List of possible mutations (ref. nucl, position, alt. nucl)
        """
        logger.debug(f"Retrieving mutations for {codon_in} -> {', '.join(codon_targets)}")
        mutations_out = []

        # Single nucleotide changes
        for idx, nucl in enumerate(codon_in):
            for nucl_alt in 'ACTG':
                if nucl == nucl_alt:
                    continue
                new_codon = codon_in[:idx] + nucl_alt + codon_in[idx + 1:]
                if new_codon not in codon_targets:
                    continue
                mutations_out.append((nucl, idx, nucl_alt))

        # Double nucleotide changes
        for i in (0, 1):
            for nucl_alt_0, nucl_alt_1 in itertools.product('ACTG', repeat=2):
                if nucl_alt_0 == codon_in[i] or nucl_alt_1 == codon_in[i + 1]:
                    continue
                codon_new = nucl_alt_0 + nucl_alt_1 + codon_in[-1] if i == 0 else codon_in[0] + nucl_alt_0 + nucl_alt_1
                if codon_new not in codon_targets:
                    continue
                mutations_out.append((codon_in[i:i + 2], i, nucl_alt_0 + nucl_alt_1))

        # Triple nucleotide changes
        for codon in codon_targets:
            mutations_out.append((codon_in, 0, codon))

        # Return mutations
        if len(mutations_out) == 0:
            logger.warning('No mutations found')
        else:
            logger.debug(f'{len(mutations_out)} mutations found')
        return mutations_out

    def __get_locus_annotation(self, gene_name: str) -> dict[str, Any]:
        """
        Returns the annotation for the input locus.
        :param gene_name: Gene name
        :return: Locus annotation
        """
        try:
            return self._gff_by_gene_name[gene_name]
        except KeyError:
            return self._gff_by_locus_tag[gene_name]

    def __parse_input_tsv_files(self) -> None:
        """
        Parses the input TSV files.
        :return: None
        """
        # Parse the TSV file with locations
        self._data_mut_locations = pd.read_table(self._args.tsv_coords)
        logger.info(f'{len(self._data_mut_locations):,} genomic coordinates parsed')

        # Parse the TSV file with the WHO AMR catalogue
        self._data_amr_catalogue = pd.read_table(self._args.tsv_catalogue)
        logger.info(f'{len(self._data_amr_catalogue):,} AMR associations parsed (WHO)')

        # Parse the TSV file with the in-house AMR catalogue
        self._data_amr_inhouse = pd.read_table(self._args.tsv_inhouse)
        logger.info(f'{len(self._data_amr_inhouse):,} AMR associations parsed (in-house)')

        # Parse antibiotics
        self._data_antibiotics = pd.read_table(self._args.tsv_abs)
        logger.info(f'{len(self._data_antibiotics):,} Antibiotics parsed')

    def __parse_input_reference(self) -> None:
        """
        Parses the input reference genome and annotation.
        :return: None
        """
        columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        data_gff = pd.read_table(self._args.ref_gff, names=columns, comment='#', header=None, sep='\t')
        logger.info(f'{len(data_gff):,} features parsed from input GFF file')

        # Retain only 'gene' features
        data_gff = data_gff[data_gff['type'] == 'gene']
        logger.info(f"{len(data_gff):,} features marked as 'gene'")

        # Store in dictionaries for retrieval
        data_gff['gene_name'] = data_gff['attributes'].apply(lambda x: re.search('gene=(.*?);', x).group(1))
        data_gff['locus_tag'] = data_gff['attributes'].apply(lambda x: re.search('locus_tag=(.*?);', x).group(1))
        self._gff_by_gene_name = {r['gene_name']: r for r in data_gff.to_dict('records')}
        self._gff_by_locus_tag = {r['locus_tag']: r for r in data_gff.to_dict('records')}

        # Parse the reference sequence
        with open(self._args.ref_fasta) as handle:
            self._seq_ref = next(SeqIO.parse(handle, 'fasta'))
        logger.info(f'Reference sequence parsed ({self._args.ref_fasta.name})')

    def __merge_amr_associations(self) -> pd.DataFrame:
        """
        Merges the in-house and WHO AMR associations.
        :return: Merged data frame
        """
        # Rename & reformat confidence column for WHO database
        self._data_amr_catalogue.rename(columns={'FINAL CONFIDENCE GRADING': 'confidence'}, inplace=True)
        self._data_amr_catalogue['confidence'] = self._data_amr_catalogue['confidence'].apply(
            lambda x: UpdateAMRDB.MAPPING_CONFIDENCE[x].value)

        # Add source column
        self._data_amr_catalogue['source'] = 'WHO'
        self._data_amr_inhouse['source'] = 'NRC'

        # Return concatenated dataframe
        target_columns = ['drug', 'gene', 'mutation', 'effect', 'confidence', 'source', 'comment']
        data_concat = pd.concat([self._data_amr_catalogue[target_columns], self._data_amr_inhouse[target_columns]])
        data_concat['variant'] = data_concat.apply(lambda row: f"{row['gene']}_{row['mutation']}", axis=1)
        return data_concat

    def __add_inhouse_positions(self) -> None:
        """
        Adds the positions of the in-house mutations to the database.
        Note that the GFF file is 1-based and the sequence is 0-based.
        :return: None
        """
        # Load the codon table and existing positions
        codon_table = CodonTable.unambiguous_dna_by_id[11]
        var_by_key = {
            (pos, ref, alt): list(d['variant']) for (pos, ref, alt), d in self._data_mut_locations.groupby(
                by=['position', 'reference_nucleotide', 'alternative_nucleotide'])}

        records_out = []
        for row in self._data_amr_inhouse.to_dict('records'):
            if row['type'] == 'PROM':
                records_out.append({
                    'gene': row['gene'],
                    'position': int(row['position']),
                    'mutation': row['mutation'],
                    'reference_nucleotide': row['reference_nucleotide'],
                    'alternative_nucleotide': row['alternative_nucleotide']
                })
            elif row['type'] == 'AA':
                # Retrieve annotation
                gene_data = self.__get_locus_annotation(row['gene'])
                logger.info(f"Processing '{row['gene']} {row['mutation']}' ({gene_data['strand']})")

                # Parse the mutation
                m = re.match(r'p.([A-z]+)(\d+)([A-z*]+)', row['mutation'])
                position = int(m.group(2))
                aa_alt = m.group(3)
                aa_alt_short = IUPACData.protein_letters_3to1.get(aa_alt)

                # Forward orientation
                if gene_data['strand'] == '+':
                    gene_start = gene_data['start']
                    codon_start = gene_start + ((position - 1) * 3)
                    ref_codon = self._seq_ref[codon_start-1:codon_start+2]
                    codons_target = [
                        str(codon) for codon, aa in codon_table.forward_table.items() if aa == aa_alt_short]
                    for n_ref, pos_rel, n_alt in UpdateAMRDB.get_mutations(str(ref_codon.seq), codons_target):
                        key_pos = (codon_start + pos_rel, n_ref, n_alt)
                        if key_pos in var_by_key:
                            if (len(var_by_key[key_pos]) > 1 or var_by_key[key_pos][0] !=
                                    f"{row['gene']}_{row['mutation']}"):
                                raise ValueError(f'Mutations at position {key_pos} do not match')
                            logger.debug(f'Skipping duplicate {key_pos}')
                            continue

                        records_out.append({
                            'gene': row['gene'],
                            'strand': gene_data['strand'],
                            'ref_codon': str(ref_codon.seq),
                            'ref_aa': codon_table.forward_table.get(str(ref_codon.seq), '-'),
                            'mutation': row['mutation'],
                            'position': codon_start + pos_rel,
                            'reference_nucleotide': n_ref,
                            'alternative_nucleotide': n_alt
                        })

                # Reverse orientation
                else:
                    gene_start = gene_data['end']
                    codon_start = gene_start - ((position - 1) * 3)
                    ref_codon = self._seq_ref[codon_start-3:codon_start]
                    ref_codon_rc = ref_codon.reverse_complement()
                    codons_target = [
                        str(codon) for codon, aa in codon_table.forward_table.items() if aa == aa_alt_short]
                    for n_ref, pos_rel, n_alt in UpdateAMRDB.get_mutations(str(ref_codon_rc.seq), codons_target):
                        records_out.append({
                            'gene': row['gene'],
                            'strand': gene_data['strand'],
                            'ref_codon': str(ref_codon.seq),
                            'ref_codon_rc': str(ref_codon_rc.seq),
                            'ref_aa': codon_table.forward_table.get(str(ref_codon_rc.seq), '-'),
                            'mutation': row['mutation'],
                            'position': codon_start - (pos_rel + len(n_alt) - 1),
                            'reference_nucleotide': str(Seq(n_ref).reverse_complement()),
                            'alternative_nucleotide': str(Seq(n_alt).reverse_complement())
                        })
            else:
                logger.warning(f"Unsupported mutation type: {row['type']}")
        data_out = pd.DataFrame(records_out)
        data_out['variant'] = data_out.apply(lambda x: f"{x['gene']}_{x['mutation']}", axis=1)
        data_out['chromosome'] = self._seq_ref.id
        target_columns = ['variant', 'chromosome', 'position', 'reference_nucleotide', 'alternative_nucleotide']
        logger.info(f'Adding {len(data_out):,} positions from in-house mutations')
        self._data_mut_locations = pd.concat([self._data_mut_locations, data_out[target_columns]])

    def __extract_gene_regions(self, data_amr_associations: pd.DataFrame) -> pd.DataFrame:
        """
        Extracts the regions with the AMR-associated genes from the parsed GFF information.
        :param data_amr_associations: Dataframe with AMR associations
        :return: Data frame with region data
        """
        records_out = []
        for gene_name, data in data_amr_associations.groupby('gene'):
            logger.debug(f"Processing '{gene_name}' ({len(data):,} associations)")

            # Get GFF information for gene
            row_gff = self.__get_locus_annotation(str(gene_name))

            # Store record
            records_out.append({
                'chr': 'Chromosome',
                'start': row_gff['start'],
                'end': row_gff['end'],
                'gene': gene_name,
                'type': re.search('gene_biotype=(.*?);', row_gff['attributes']).group(1),
                'abs': ', '.join(sorted(data['drug'].unique()))
            })
        return pd.DataFrame(records_out)

    def __extract_prom_regions(self, data_amr_associations: pd.DataFrame) -> pd.DataFrame:
        """
        Extracts the regions with the promoters of the AMR-associated genes that contain mutations from the catalogue.
        :param data_amr_associations: Dataframe with AMR associations
        :return: Data frame with region data
        """
        # Merge data frames (note that a single variant can match multiple positions)
        data_merged = pd.merge(data_amr_associations, self._data_mut_locations, how='left', on='variant')

        # Identify promotor mutations
        data_merged['is_promotor'] = data_merged['mutation'].apply(lambda x: x.startswith('c.-') or x.startswith('n.-'))

        # Add regions
        records_out = []
        for gene, data_mutations in data_merged[data_merged['is_promotor']].groupby('gene'):
            p_start, p_end = (int(data_mutations['position'].min()), int(data_mutations['position'].max()))
            logger.debug(f'{len(data_mutations):,} promotor mutations for {gene}')
            records_out.append({
                'chr': 'Chromosome',
                'start': p_start,
                'end': p_end,
                'gene': f'{gene}_prom',
                'type': 'intergenic',
                'abs': ', '.join(sorted(data_mutations['drug'].unique()))
            })
        return pd.DataFrame(records_out)

    @staticmethod
    def __get_matching_region(position: int, regions: list[dict]) -> Union[str, None]:
        """
        Returns the region matching that covers the input position (if available).
        :param position: Genome position
        :param regions: Extracted genomic regions
        :return: The matching region (if available), None otherwise
        """
        for row in regions:
            if row['start'] <= position <= row['end']:
                return row['gene']
        return None

    def __verify_db_content(self, data_regions_all: pd.DataFrame, data_amr_associations: pd.DataFrame) -> None:
        """
        Verifies the database content and checks which mutations are covered.
        :param data_regions_all: Dataframe containing AMR regions
        :param data_amr_associations: Dataframe containing AMR associations
        :return: None
        """
        # Extract mutations that fall outside the regions in the BED file
        regions_list = data_regions_all.to_dict('records')
        self._data_mut_locations['region'] = self._data_mut_locations['position'].apply(
            lambda x: UpdateAMRDB.__get_matching_region(x, regions_list))
        mutations_not_covered = self._data_mut_locations[pd.isna(self._data_mut_locations['region'])]
        logger.info(f"{len(mutations_not_covered):,}/{len(self._data_mut_locations):,} mutations not covered")
        mutations_not_covered.to_csv(self._args.dir_out / '_not_covered.tsv', sep='\t', index=False)

        # Check the possible mutation positions for the AMR associations
        muts_by_variant = {variant: data for variant, data in self._data_mut_locations.groupby('variant')}
        data_amr_associations['nb_muts'] = data_amr_associations['variant'].apply(
            lambda x: len(muts_by_variant.get(x, [])))
        logger.info(f"No location found for {sum(data_amr_associations['nb_muts'] == 0)} mutations")
        with open(self._args.dir_out / '_not_in_locations.tsv', 'w') as handle_out:
            for variant in sorted(list(data_amr_associations[data_amr_associations['nb_muts'] == 0]['variant'])):
                handle_out.write(variant)
                handle_out.write('\n')

        # Verify that the antibiotics in the AMR catalogue are present in the TSV file
        short_by_full = {ab['AB']: ab for ab in self._data_antibiotics.to_dict('records')}
        data_amr_associations['drug'].apply(lambda x: short_by_full[x])
        logger.info('All antibiotics from AMR catalogue are present in the input TSV file')

    def __export_db_files(self, data_amr_associations: pd.DataFrame) -> None:
        """
        Exports the database files to the output directory.
        :param data_amr_associations: Dataframe containing AMR associations
        :return: None
        """
        # Copy the TSV files to the output folder
        self._data_mut_locations.to_csv(self._args.dir_out / 'mutation_locations.tsv', sep='\t', index=False)
        data_amr_associations.to_csv(self._args.dir_out / 'amr_associations_all.tsv', sep='\t', index=False)
        shutil.copyfile(self._args.tsv_abs, self._args.dir_out / 'antibiotics.tsv')

        # Create a text file with the DB version
        with open(self._args.dir_out / 'VERSION', 'w') as handle:
            handle.write(self._args.version)
            handle.write('\n')
        logger.info('Version file created')

    def run(self) -> None:
        """
        Runs the updater.
        :return: None
        """
        # Parse input files
        self.__parse_input_tsv_files()
        self.__parse_input_reference()

        # Merge the WHO and in-house mutations
        data_amr_associations = self.__merge_amr_associations()
        self.__add_inhouse_positions()

        # Construct BED file with target gene regions
        data_regions_genes = self.__extract_gene_regions(data_amr_associations)
        logger.info(f'{len(data_regions_genes):,} gene regions extracted')

        # As the promoter mutations fall outside the gene regions, they need to be extracted separately
        data_regions_prom = self.__extract_prom_regions(data_amr_associations)
        logger.info(f'{len(data_regions_genes):,} promotor regions extracted')

        # Export merged BED file
        data_regions_all = pd.concat([data_regions_genes, data_regions_prom])
        data_regions_all.sort_values(by='start', inplace=True)
        data_regions_all.to_csv(self._args.dir_out / 'amr_regions.bed', sep='\t', header=False, index=False)
        logger.info('BED file created')

        # Verify & export output
        self.__verify_db_content(data_regions_all, data_amr_associations)
        self.__export_db_files(data_amr_associations)


if __name__ == '__main__':
    initialize_logging()
    updater = UpdateAMRDB()
    updater.run()
