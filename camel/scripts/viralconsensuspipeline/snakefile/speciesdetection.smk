from pathlib import Path

from camel.app.core.snakemake import snakemakeutils


rule species_detection_report:
    """
    Creates the HTML output report for the species detection.
    """
    input:
        TSV = config['species_detection']['tsv'] if config.get('species_detection') is not None else []
    output:
        HTML = 'species_detection/report/html.iob' # speciesdetection.OUTPUT_REPORT
    params:
        group = (config.get('species_detection') or {}).get('group'),
        species = (config.get('species_detection') or {}).get('species')
    run:
        import pandas as pd
        from camelcore.app.io.tooliovalue import ToolIOValue
        from camelcore.app.reports.htmlreportsection import HtmlReportSection

        data_species = pd.read_table(input.TSV)
        section = HtmlReportSection('Species detection')
        section.add_paragraph(
            f'The species was determined by competitive read mapping against the reference genomes for: '
            f'<i>{params.group}</i>. Detected species: <b>{params.species}</b>')
        section.add_table(
            [[row['species'], f"{row['mapped_bases']:,}", f"{100 * row['fraction_bases']:.2f}",
              f"{100 * row['covered_rate']:.2f}"] for row in data_species.to_dict('records')],
            ['Species', 'Mapped bases', 'Mapped bases (%)', 'Covered (%)'],
            [('class', 'data')])
        snakemakeutils.dump_object([ToolIOValue(section)], Path(output.HTML))

rule species_detection_informs:
    """
    Collects the informs of the species detection for the commands section.
    """
    input:
        JSON = config['species_detection']['informs'] if config.get('species_detection') is not None else []
    output:
        INFORMS = 'species_detection/informs/informs.io' # speciesdetection.OUTPUT_INFORMS
    run:
        import json
        with open(input.JSON) as handle:
            informs_all = json.load(handle)
        for inform in informs_all:
            inform['_tag'] = 'Species detection'
        snakemakeutils.dump_object(informs_all, Path(output.INFORMS))

rule species_detection_summary:
    """
    Creates the summary output for the species detection.
    """
    input:
        TSV = config['species_detection']['tsv'] if config.get('species_detection') is not None else []
    output:
        FILE = 'species_detection/summary/summary.{ext}' # speciesdetection.OUTPUT_SUMMARY
    params:
        species = (config.get('species_detection') or {}).get('species'),
        ext = lambda wildcards: wildcards.ext
    run:
        import pandas as pd
        data_species = pd.read_table(input.TSV)
        data_summary = [('species_detection_species', params.species)]
        for row in data_species.to_dict('records'):
            data_summary.append((f"species_detection_{row['species']}_fraction_bases", row['fraction_bases']))
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'species_detection')
