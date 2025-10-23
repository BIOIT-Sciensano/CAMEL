import json
from pathlib import Path

from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import assembly
from camel.scripts.salmonellapipeline.snakefile import serotyping_sistr


rule serotyping_sistr_run:
    """
    Runs SISTR in cgMLST mode.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA
    output:
        JSON = 'serotyping/sistr/tool/json.io', # serotyping_sistr.OUTPUT_JSON
        INFORMS = 'serotyping/sistr/tool/informs.io' # serotyping_sistr.OUTPUT_INFORMS
    params:
        dir_ = 'serotyping/sistr/tool',
        db = config['serotyping']['sistr']['path']
    run:
        from camel.app.tools.pipelines.salmonella.sistr import Sistr
        sistr_tool = Sistr()
        sistr_tool.add_input_files({'DIR': [ToolIODirectory(Path(params.db))]})
        snakemakeutils.add_pickle_inputs(sistr_tool, input)
        step = Step(rule_name=str(rule), tool=sistr_tool, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(sistr_tool, output)

rule serotyping_sistr_dump_summary_info:
    """
    Creates the summary output for SISTR.
    """
    input:
        JSON = rules.serotyping_sistr_run.output.JSON,
        INFORMS = rules.serotyping_sistr_run.output.INFORMS
    output:
        FILE = 'serotyping/sistr/summary/summary_out_sistr.{ext}' # serotyping_sistr.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        # Parse the SISTR output
        with snakemakeutils.load_object(Path(input.JSON))[0].path.open('r') as handle:
            json_data = json.load(handle)[0]

        # Reformat data
        header_locus = ['Locus', 'serotype_or_group', '% Identity', 'HSP/Locus length', 'Contig', 'Position in contig']
        if json_data['qc_status'] == 'PASS':
            hits_dict_tsv = {
                'serotype_antigenic_formula':':'.join([
                    str(json_data['o_antigen']), str(json_data['h1']), str(json_data['h2'])]),
                'serotype_serogroup': json_data['serogroup'],
                'serotype_consensus': json_data['serovar'],
                'qc_status' : 'PASS'
            }
            serotyping_sistr.sistr_output_parser(json_data['h1_flic_prediction'], 'fliC', 'h1', hits_dict_tsv)
            serotyping_sistr.sistr_output_parser(json_data['h2_fljb_prediction'], 'fljB', 'h2', hits_dict_tsv)
            serotyping_sistr.sistr_output_parser(json_data['serogroup_prediction']['wzx_prediction'], 'wzx', 'o', hits_dict_tsv)
            serotyping_sistr.sistr_output_parser(json_data['serogroup_prediction']['wzy_prediction'], 'wzy', 'o', hits_dict_tsv)
        else:
            hits_dict_tsv = {
                'serotype_antigenic_formula': '-',
                'serotype_serogroup': '-',
                'serotype_consensus': '-',
                'qc_status': 'FAIL'}
            for variable in ['hits_serotype_h1_fliC', 'hits_serotype_h2_fljB', 'hits_serotype_o_wzx', 'hits_serotype_o_wzy']:
                hits_dict_tsv[variable] = '-'

        # Tool information
        informs_sistr = snakemakeutils.load_object(Path(str(input.INFORMS)))
        rows_out = [(f'sistr_{key}', value) for key, value in hits_dict_tsv.items()]
        rows_out.extend([
            (f'sistr_tool_version', informs_sistr['_name_full']),
            (f'sistr_db_version', informs_sistr['last_update_date'])
        ])

        # Create JSON output
        if params.ext == 'json':
            entries = ['sistr_hits_serotype_h1_fliC', 'sistr_hits_serotype_h2_fljB', 'sistr_hits_serotype_o_wzx', 'sistr_hits_serotype_o_wzy']
            rows_out = [(k, snakemakeutils.convert_list_to_dict([v.split(',')], header_locus)) if k in entries and v != '-' else (k, v) for k, v in rows_out]
        snakemakeutils.export_summary(rows_out, Path(output.FILE), str(params.ext), 'sistr')

rule serotyping_sistr_report:
    """
    This rule creates a simple output report for the SISTR serotyping analysis.
    """
    input:
        JSON_SISTR = rules.serotyping_sistr_run.output.JSON,
        INFORMS_serotyping_sistr = rules.serotyping_sistr_run.output.INFORMS
    output:
        VAL_HTML = 'serotyping/sistr/report/html.iob' # serotyping_sistr.OUTPUT_REPORT
    params:
        dir_ = 'serotyping/sistr/report' ,
        db_path_sistr = config['serotyping']['sistr']['path']
    run:
        from camel.app.tools.pipelines.salmonella.sistrreporter import SistrReporter
        reporter = SistrReporter()
        reporter.add_input_files({'DIR_sistr': [ToolIODirectory(Path(params.db_path_sistr))]})
        snakemakeutils.add_pickle_inputs(reporter, input)
        step = Step(str(rule), reporter, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule serotyping_sistr_report_empty:
    """
    Creates an empty HTML report for the SISTR serotyping analysis.
    """
    output:
        VAL_HTML = 'serotyping/sistr/report/html-empty.iob' # serotyping_sistr.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.tools.pipelines.salmonella.sistrreporter import SistrReporter
        snakepipelineutils.create_empty_report_section(SistrReporter.TITLE, Path(output.VAL_HTML))
