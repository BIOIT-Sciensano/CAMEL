import json
from pathlib import Path

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class SeqSero2Reporter(Tool):
    """
    Parses Seqsero2's TSV output results and returns an HTML report.
    """
    TITLE = 'SeqSero2'

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('SeqSero2 Reporter', '0.1')
        self._section = None

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._section = HtmlReportSection(
            SeqSero2Reporter.TITLE,subtitle=self._input_informs['serotyping_seqsero2']['_name_full'])
        self.__add_section_seqsero()
        self.__add_file_output()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_section_seqsero(self) -> None:
        """
        Adds the SeqSero2 section to the HTML report.
        :return: None
        """
        self._section.add_header(self._input_informs['serotyping_seqsero2']['_name'], 2)

        # Mandatory seqsero assay (from FASTA)
        self._section.add_header('SeqSero2 serotyping - assembly kmer mode', 4)
        self.___add_table_serotype_seqsero(self._tool_inputs['TXT_seqsero2_kmer'][0].path)

        # Optional seqsero assays (from FASTQ)
        self._section.add_header('SeqSero2 serotyping - raw read allele mode', 4)
        if 'TXT_seqsero2_allele' in self._tool_inputs:
            self.___add_table_serotype_seqsero(self._tool_inputs['TXT_seqsero2_allele'][0].path)
        else:
            self._section.add_paragraph('SeqSero2 serotyping (raw read allele mode) not available in FASTA-input mode')
        self._section.add_header('SeqSero2 serotyping - raw read kmer mode', 4)
        if 'TXT_seqsero2_kmerread' in self._tool_inputs:
            self.___add_table_serotype_seqsero(self._tool_inputs['TXT_seqsero2_kmerread'][0].path)
        else:
            self._section.add_paragraph('SeqSero2 serotyping (assembly kmer mode) not available in FASTA-input mode')

        # add last update of the seqsero2 db
        db_dir = self._tool_inputs['DIR_seqsero2'][0].path
        self.___add_database_information(db_dir)

    def ___add_table_serotype_seqsero(self, input_file_path: Path) -> None:
        """
        Generates and adds the table for seqsero2 tool.
        :param input_file_path: the text file containing the results for a seqsero2 run in any mode.
        :return: None
        """
        resultsdict = {}
        with input_file_path.open('r') as handle:
            for line in handle:
                parts = line.rstrip().split('\t')
                resultsdict[parts[0]] = parts[1] if len(parts) > 1  else ''
        table_data = []
        header = ['O-antigen', 'H1-antigen (fliC)', 'H2-antigen (fljB)', 'Antigenic formula', 'Serotype']
        row = [resultsdict['O antigen prediction:'], resultsdict['H1 antigen prediction(fliC):'],
               resultsdict['H2 antigen prediction(fljB):'], resultsdict['Predicted antigenic profile:'],
               resultsdict['Predicted serotype:']]
        table_data.append(row)
        self._section.add_table(table_data, header, [('class', 'data')])

    def __add_file_output(self) -> None:
        """
        Add the output tsv file to the html.
        :return: None
        """
        # Mandatory seqsero assay (from FASTA)
        relative_path = Path('serotyping', 'seqsero2', 'summary_out_seqsero2_kmer.tsv')
        self._section.add_link_to_file("Download SeqSero2 Kmer mode (TSV)", relative_path)
        self._section.add_file(self._tool_inputs['TXT_seqsero2_kmer'][0].path, relative_path)

        # Optional seqsero assays (from FASTQ)
        if 'TXT_seqsero2_allele' in self._tool_inputs:
            relative_path = Path('serotyping', 'seqsero2', 'summary_out_seqsero2_allele.tsv')
            self._section.add_link_to_file("Download SeqSero2 Allele mode (TSV)", relative_path)
            self._section.add_file(self._tool_inputs['TXT_seqsero2_allele'][0].path, relative_path)
        if 'TXT_seqsero2_kmerread' in self._tool_inputs:
            relative_path = Path('serotyping', 'seqsero2', 'summary_out_seqsero2_kmerread.tsv')
            self._section.add_link_to_file("Download SeqSero2 Kmer reads mode (TSV)", relative_path)
            self._section.add_file(self._tool_inputs['TXT_seqsero2_kmerread'][0].path, relative_path)

    def ___add_database_information(self, db_dir: Path) -> None:
        """
        Adds the date of latest database update.
        :param db_dir: Input database directory
        :return: None
        """
        db_metadata_file = db_dir / 'db_update_info.json'
        if not db_metadata_file.is_file():
            raise FileNotFoundError(f'Database metadata not found: {db_metadata_file}')
        with db_metadata_file.open() as handle:
            metadata = json.load(handle)
            last_update_date = metadata['last_update_date']
        self._section.add_paragraph(f'Last updated: {last_update_date}')
