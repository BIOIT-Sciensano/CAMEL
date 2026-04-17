from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils


rule seqtk_seq:
    input:
        FASTA = config['input']
    output:
        FASTA = 'seqtk_seq/fasta.io',
        INFORMS = 'seqtk_seq/informs.io'
    params:
        dir_ = 'seqtk_seq'
    run:
        from camel.app.tools.seqtk.seqtkseq import SeqtkSeq
        seqtk_seq_tool = SeqtkSeq()
        seqtk_seq_tool.add_input_files({
            'FASTA': [ToolIOFile(Path(input.FASTA))]
        })
        seqtk_seq_tool.update_parameters(output_filename='filtered.fasta')
        step = Step(rule_name=str(rule), tool=seqtk_seq_tool, dir_=Path(params.dir_).absolute())
        step.run()
        snakemakeutils.dump_io_outputs(seqtk_seq_tool, output)

rule seqkit_stats:
    input:
        FASTA = rules.seqtk_seq.output.FASTA
    output:
        INFORMS = 'seqkit_stats/informs.io'
    params:
        dir_ = 'seqkit_stats'
    run:
        from camel.app.tools.seqkit.seqkitstats import SeqkitStats
        seqtk_stats = SeqkitStats()
        snakemakeutils.add_io_inputs(seqtk_stats, input)
        step = Step(rule_name=str(rule), tool=seqtk_stats, dir_=Path(params.dir_).absolute())
        step.run()
        snakemakeutils.dump_io_outputs(seqtk_stats, output)

rule collect_output:
    """
    Collects the tool output.
    """
    input:
        INFORMS = rules.seqkit_stats.output.INFORMS
    output:
        TSV = config['output']
    run:
        informs = snakemakeutils.load_object(Path(input.INFORMS))
        with open(output.TSV, 'w') as handle:
            handle.write('\t'.join(['total_length', str(informs['sum_len'])]))
            handle.write('\n')
            handle.write('\t'.join(['version', informs['_version']]))
            handle.write('\n')
