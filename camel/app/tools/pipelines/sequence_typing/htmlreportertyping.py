import json
from pathlib import Path
from typing import Optional

from Bio import SeqIO

from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlexpandablediv import HtmlExpandableDiv
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.components.sequencetyping.typinghitbase import TypingHitBase
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


class HtmlReporterTyping(Tool):
    """
    Tool that creates HTML reports for the sequence typing pipeline.

    Input:
        - HTML: Path to the HTML file to write the report
        - DIR: Directory to store files that are included in the HTML report
        - Informs 'Scheme': Information about the scheme
        - VAL_Hits: Hits detected for each locus
    Output:
        - HTML: Path to the generated report
    """

    INFO_FILENAME = 'sequence_typing.json'
    NEW_ALLELES_FILENAME = 'new_alleles.fasta'

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('HTML Reporter', '0.1')
        self._report_section = None
        self._output_folder = None
        self._sub_folder = None

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        import pprint
        pprint.pprint(self._input_informs)

        self.__initialize_report()
        if 'ST' in self._input_informs:
            self.__add_sequence_type()

        # Add tables with the detected hits
        add_subtitle = True if all(len(self._tool_inputs[f'hits_{key}']) > 1 for key in ('nucl', 'pept')) else False
        if len(self._tool_inputs['hits_nucl']) != 0:
            self.__add_output_table(self._tool_inputs['TSV_nucl'][0].path, self._tool_inputs['hits_nucl'],
                                    'Nucleotide loci' if add_subtitle else None)
        if len(self._tool_inputs['hits_pept']) != 0:
            self.__add_output_table(self._tool_inputs['TSV_pept'][0].path, self._tool_inputs['hits_pept'],
                                    'Peptide loci' if add_subtitle else None)

        if 'forced_detection_method' in self._parameters:
            self._report_section.add_alert(
                f"Allele detection performed with <b>{self._parameters['forced_detection_method'].value}</b>.", 'info')
        self._add_novel_alleles_section()
        self.add_scheme_info_section()

        # Add a custom message (if specified in the parameters)
        if 'message' in self._parameters:
            self._report_section.add_alert(
                self._parameters['message'].value, self._parameters['message_category'].value)
        
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section)]
        self.__export_analysis_metadata()

    def __initialize_report(self) -> None:
        """
        Initializes the HTML report.
        :return: None
        """
        self._report_section = HtmlReportSection(self._input_informs['scheme']['title'], 3)
        self._sub_folder = Path('sequence_typing', FileSystemHelper.make_valid(self._input_informs['scheme']['name']))

    def __add_sequence_type(self) -> None:
        """
        Adds the sequence type to the report.
        :return: None
        """
        profile_data = self._input_informs['ST']
        header = ['Loci matched', '% matched'] + [key for key, _ in profile_data['metadata']]
        # Determine color for ST cell
        if profile_data['is_detected'] and profile_data['percent_detected'] == 100:
            color = 'green'
        elif profile_data['is_detected']:
            color = 'yellow'
        else:
            color = 'red'
        table_data = [
            [f"{profile_data['nb_detected']}/{profile_data['nb_loci']}",
             f"{profile_data['percent_detected']:.2f}%",
             HtmlTableCell(profile_data['symbol'], color)] +
            [value if value != '' else '-' for _, value in profile_data['metadata'][1:]]]

        self._report_section.add_table(table_data, header, table_attributes=[('class', 'data')])

    def __add_output_table(self, output_tsv: Path, hits_io: list[ToolIOValue], sub_header: Optional[str]) -> None:
        """
        Adds the output table with the detected alleles.
        :param output_tsv: Tabular output file
        :param hits_io: Detected hits
        :param sub_header: If not None, this sub header is added to the report
        :return: None
        """
        table_header = hits_io[0].value.html_column_names()
        table_data = [h.value.to_html_row(self._report_section, self._sub_folder) for h in sorted(
            hits_io, key=lambda x: x.value.locus)]

        if sub_header is not None:
            self._report_section.add_header(sub_header, 4)

        # Add slider for big tables
        is_hidden = len(hits_io) > 12 or ('hidden' in self._parameters)
        if is_hidden:
            div = HtmlExpandableDiv('table-{}'.format(
                self._input_informs['scheme']['name'].lower()), f'alleles ({len(hits_io)})')
            div.add_table(table_data, table_header, [('class', 'data')])
            self._report_section.add_html_object(div)
        else:
            self._report_section.add_table(table_data, table_header, [('class', 'data')])

        # Add regular TSV file
        relative_path = self._sub_folder / output_tsv.name
        self._report_section.add_file(output_tsv, relative_path)
        self._report_section.add_link_to_file("Download (TSV)", relative_path)

        # Add TSV with hashes (if it exists)
        output_tsv_hashes = output_tsv.parent / output_tsv.name.replace('.tsv', '-hashes.tsv')
        if not output_tsv_hashes.exists():
            return
        relative_path_hashes = self._sub_folder / output_tsv_hashes.name
        self._report_section.add_file(output_tsv_hashes, relative_path_hashes)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'scheme' not in self._input_informs:
            raise InvalidToolInputError("Scheme information is required")
        super()._check_input()

    def __export_analysis_metadata(self) -> None:
        """
        Exports the analysis metadata file. The information can be used for further processing of the sequence typing
        output (e.g., for generating MLST trees).
        :return: None
        """
        path = self.folder / HtmlReporterTyping.INFO_FILENAME
        with path.open('w') as handle:
            json.dump({
                'scheme': self._input_informs['scheme']['name'],
                'sample': self._tool_inputs['VAL_SAMPLE'][0].value
            }, handle)
        self._report_section.add_file(path, self._sub_folder / HtmlReporterTyping.INFO_FILENAME)

    def add_scheme_info_section(self) -> None:
        """
        Adds the section with information about the typing scheme.
        :return: None
        """
        self._report_section.add_header('Scheme info', level=4)
        table_data = [
            ('Last scheme update', self._input_informs['scheme']['update']['last_update_date']),
            ('Last scheme change', self._input_informs['scheme']['update'].get('last_change_date')),
            ('Origin', self._input_informs['scheme']['update'].get('origin')),
        ]
        for i in range(0, len(table_data)):
            if table_data[i][-1] is not None:
                continue
            table_data[i] = (table_data[i][0], 'n/a')
        self._report_section.add_table(table_data, ['Field', 'Value'], [('class', 'data')])

    def __generate_novel_allele_fasta(self, hits_novel: list[TypingHitBase]) -> Path:
        """
        Generates a FASTA file with the novel alleles.
        :param hits_novel: Hits with novel alleles
        :return: Path to FASTA file
        """
        # Collect sequences
        seqs_out = []
        for hit in hits_novel:
            seqs_out.append(SeqIO.SeqRecord(
                seq=hit.new_allele_sequence,
                id=hit.locus,
                description='potential novel allele'
            ))

        # Save in FASTA format
        basename = f"novel-{self._tool_inputs['VAL_SAMPLE'][0].value}-{self._input_informs['scheme']['name']}.fasta"
        path_out = self.folder / FileSystemHelper.make_valid(basename)
        with path_out.open('w') as handle:
            SeqIO.write(seqs_out, handle, 'fasta')
        return path_out

    def _add_novel_alleles_section(self) -> None:
        """
        Adds the section with the novel alleles (if there are any).
        :return: None
        """
        # Check if there are (potential) novel alleles
        hits_all = [h.value for h in self._tool_inputs.get('hits_nucl', []) + self._tool_inputs.get('hits_pept', [])]
        hits_novel_allele = [h for h in hits_all if h.is_new_allele()]
        if len(hits_novel_allele) == 0:
            return
        logger.info(f'{len(hits_novel_allele)} potential novel alleles detected')

        # Create a FASTA file and add it to the report
        path_fasta = self.__generate_novel_allele_fasta(hits_novel_allele)
        relative_path = self._sub_folder / path_fasta.name

        # Add the novel allele information
        self._report_section.add_header('Novel alleles (*)', level=4)
        self._report_section.add_paragraph('Potentially novel alleles were detected for the following loci: {}'.format(
            ', '.join(f'<b>{h.locus}</b>' for h in hits_novel_allele)))
        self._report_section.add_paragraph(
            'These allele sequences can be submitted to the corresponding database after verification (intact open-'
            'reading frames, etc). Note that not all novel alleles might have been detected by this screening.')
        self._report_section.add_file(path_fasta, relative_path)
        self._report_section.add_link_to_file('Download consensus sequence(s) (FASTA)', relative_path)
