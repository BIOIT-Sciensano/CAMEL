import logging
from pathlib import Path

import pandas as pd

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils


rule multi_allelic_sites_bcftools_pileup:
    """
    Creates a pileup using BCFtools.
    """
    input:
        FASTA = 'iterative_mapping/output/fasta.io',
        BAM = 'iterative_mapping/output/bam.io'
    output:
        VCF = 'multi_allelic/pileup/vcf.io'
    params:
        dir_ = 'multi_allelic/pileup',
        max_depth = 25000,
        config = 'illumina' if config['input']['type'] == 'illumina' else 'ont'
    run:
        from camel.app.tools.bcftools.bcftoolsmpileup import BcftoolsMpileup
        mpileup = BcftoolsMpileup()
        snakemakeutils.add_pickle_inputs(mpileup, input)
        mpileup.update_parameters(
            max_depth=params.max_depth, output_type='v', annotate='FORMAT/AD', config=params.config)
        step = Step(rule_name=str(rule), tool=mpileup, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(mpileup, output)

rule multi_allelic_sites_extract_snps:
    """
    Extracts SNPs from the input pileup. 
    """
    input:
        VCF = rules.multi_allelic_sites_bcftools_pileup.output.VCF
    output:
        VCF = 'multi_allelic/snps/vcf.io'
    params:
        dir_ = 'multi_allelic/snps'
    run:
        from camel.app.tools.bcftools.bcftoolsview import BcftoolsView
        bcftools_view = BcftoolsView()
        snakemakeutils.add_pickle_inputs(bcftools_view, input)
        bcftools_view.update_parameters(types='snps')
        step = Step(rule_name=str(rule), tool=bcftools_view, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(bcftools_view, output)

rule multi_allelic_sites_extract_sites:
    """
    Extracts the multi-allelic sites from the input file with SNPs.
    """
    input:
        VCF = rules.multi_allelic_sites_extract_snps.output.VCF
    output:
        TSV = 'multi_allelic/extract_sites/tsv.io',
        INFORMS = 'multi_allelic/extract_sites/informs.io'
    params:
        dir_ = 'multi_allelic/extract_sites',
        min_dp = 20,
        min_freq_minor_allele = 0.4
    run:
        from camel.app.tools.pipelines.viral_consensus.callmultiallelicsites import CallMultiAllelicSites
        call_multi_allelic = CallMultiAllelicSites()
        snakemakeutils.add_pickle_inputs(call_multi_allelic, input)
        step = Step(rule_name=str(rule), tool=call_multi_allelic, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(call_multi_allelic, output)

rule multi_allelic_sites_apply_to_consensus:
    """
    Updates the consensus sequence by replacing the nucleotides at the multi-allelic sites with the corresponding IUPAC 
    codes.
    """
    input:
        TSV = rules.multi_allelic_sites_extract_sites.output.TSV,
        FASTA = 'iterative_mapping/output/fasta.io'
    output:
        FASTA = 'multi_allelic/updated_consensus/fasta.io'
    params:
        dir_ = 'multi_allelic/updated_consensus'
    run:
        from Bio import SeqIO
        from Bio.Seq import Seq
        from camel.app.core.io.tooliofile import ToolIOFile

        # Parse the input FASTA file
        path_fasta = snakemakeutils.load_object(Path(input.FASTA))[0].path
        with path_fasta.open() as handle:
            seq_by_id = {s.id: str(s.seq).upper() for s in SeqIO.parse(handle, 'fasta')}

        # Parse the multi-allelic sites
        path_tsv = snakemakeutils.load_object(Path(input.TSV))[0].path

        try:
            data_multi = pd.read_table(path_tsv)
            multi_allelic_sites_by_seq_id = {
                seq_id: data.to_dict('records') for seq_id, data in data_multi.groupby('chrom')}
        except pd.errors.EmptyDataError:
            logging.warning('No multi-allelic sites found')
            multi_allelic_sites_by_seq_id = {}

        # Apply the variants
        seq_new_by_id = {id_: seq_in for id_, seq_in in seq_by_id.items()}
        for id_, seq_in in seq_by_id.items():
            seq_new = seq_in
            # Apply ambiguous bases
            if id_ not in multi_allelic_sites_by_seq_id:
                continue
            for row in multi_allelic_sites_by_seq_id[id_]:
                seq_new = seq_new[:row['pos']-1] + row['iupac'] + seq_new[row['pos']:]
            seq_new_by_id[id_] = seq_new

        # Save sequences
        path_out = Path(params.dir_) / f"{path_fasta.name.replace('.fasta', '')}-ambiguous.fasta"
        with path_out.open('w') as handle:
            seqs_out = [SeqIO.SeqRecord(Seq(seq), id=id_, description='') for id_, seq in sorted(seq_new_by_id.items())]
            SeqIO.write(seqs_out, handle, 'fasta')
        snakemakeutils.dump_object([ToolIOFile(path_out)], Path(output.FASTA))

rule multi_allelic_sites_report:
    """
    Creates the output report with the multi-allelic sites.
    """
    input:
        FASTA = rules.multi_allelic_sites_apply_to_consensus.output.FASTA,
        TSV = rules.multi_allelic_sites_extract_sites.output.TSV,
        INFORMS_calling = rules.multi_allelic_sites_extract_sites.output.INFORMS
    output:
        VAL_HTML = 'multi_allelic/report/html.iob' # multiallelic.OUTPUT_REPORT
    params:
        dir_ = 'multi_allelic/report',
        name = config['input']['sample_name']
    run:
        from camel.app.tools.pipelines.viral_consensus.reportmultiallelic import ReporterMultiAllelic
        reporter = ReporterMultiAllelic()
        snakemakeutils.add_pickle_inputs(reporter, input)
        reporter.update_parameters(name=params.name)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule multi_allelic_create_summary:
    """
    Creates the summary output for the multi-allelic site calling.
    """
    input:
        INFORMS = rules.multi_allelic_sites_extract_sites.output.INFORMS
    output:
        TSV = 'multi_allelic/report/summary.{ext}' # multiallelic.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs = snakemakeutils.load_object(Path(input.INFORMS))
        data_summary = [('multi_allelic_sites_nb', int(informs['nb_sites']))]
        snakemakeutils.export_summary(data_summary, Path(output.TSV), str(params.ext), 'multiallelic')
