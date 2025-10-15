import datetime
from pathlib import Path
from typing import Any

import pandas as pd

from camel.app.camel import Camel
from camel.app.components.html import PATH_JQUERY
from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.config import config
from camel.app.error import InvalidToolInputError
from camel.app.tools.tool import Tool
from camel.resources import CSS_STYLE


class HybridAssemblyReporter(Tool):
    """
    This tool is used to generate HTML report sections based on the hybrid assembly pipeline output.
    """

    TITLE = 'Hybrid assembly pipeline 0.1'
    MATCH_COLORS = {0: None, 1: 'grey', 2: 'lightgreen', 3: 'green'}
    REPORT_STRUCTURE = [
        ['Read trimming and basic QC', 'trim'],
        ['Assembly statistics', 'assembly'],
        ['Variant calling analysis', 'vc'],
        ['Commands', 'commands'],
        ['Citations', 'citations']
    ]

    def __init__(self) -> None:
        """
        Initializes the tool.
        :return: None
        """
        super().__init__('Hybrid assembly pipeline reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid.
        :return: None
        """
        # Check tool input files
        required_inputs = [
            'FASTA', 'TSV_quast', 'TSV_vc', 'TSV_mapping', 'TSV_sniffles', 'VCF_sniffles', 'HTML_trim_illumina',
            'HTML_trim_ont', 'HTML_quast', 'WIGGLE_ale']
        for key in required_inputs:
            if key not in self._tool_inputs:
                raise InvalidToolInputError(f"Required tool input '{key}' is missing")

        # Check informs
        required_informs = [
            'quast', 'freebayes', 'clair3', 'mapping', 'sniffles', 'ale', 'sample_name', 'pipeline', 'input',
            'citations']
        for key in required_informs:
            if key not in self._input_informs:
                raise InvalidToolInputError(f"Required inform '{key}' is missing")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        path_out = Path(self._parameters['output_filename'].value)
        self._output_dir = Path(self._parameters['output_dir'].value)

        # Initialize report
        self.report = HtmlReport(path_out, self._output_dir, [PATH_JQUERY])
        self.report.initialize('Hybrid assembly pipeline - 0.1', CSS_STYLE)
        self.report.add_pipeline_header(HybridAssemblyReporter.TITLE)

        # Add content sections
        self.__add_input_section()
        self.__add_overview_links()
        self.report.add_module_header('Overview')
        self.__add_overview_section()

        # Trimming section
        self.report.add_module_header('Read trimming', id_='trim')
        self.report.add_html_object(self._tool_inputs['HTML_trim_illumina'][0].value)
        self._tool_inputs['HTML_trim_illumina'][0].value.copy_files(self.report.output_dir)
        self.report.add_html_object(self._tool_inputs['HTML_trim_ont'][0].value)
        self._tool_inputs['HTML_trim_ont'][0].value.copy_files(self.report.output_dir)

        # Other section
        self.report.add_module_header('Assembly statistics', id_='assembly')
        self.__add_quast_section(self._tool_inputs['TSV_quast'][0].path, self._input_informs['quast'])
        self.report.add_module_header('Variant calling', id_='vc')
        self.__add_vc_table(self._tool_inputs['TSV_vc'][0].path)
        self.__add_sniffles_table(self._tool_inputs['TSV_sniffles'][0].path)
        self.__add_mapping_table(self._tool_inputs['TSV_mapping'][0].path)
        self.__add_ale_score_table()

        # Commands & citations
        self.report.add_module_header('Commands', id_='commands')
        self.report.add_html_object(self._input_informs['commands'][0].value)
        self.report.add_module_header('Citations', id_='citations')
        self.report.add_html_object(self._input_informs['citations'][0].value)
        self.report.save()

    def __add_input_section(self) -> None:
        """
        Adds the information about the input data.
        :return: None
        """
        input_files = ReportPipeline.format_input_string(self._input_informs['input'])
        table_data = [
            ['Sample:', self._input_informs['sample_name']],
            ['Analysis date:', datetime.datetime.now().strftime(config.date_fmt)],
            ['Pipeline version:', self._input_informs['pipeline']['version']],
            ['Input files:', input_files]
        ]
        section = HtmlReportSection('Input')
        section.add_table(table_data, table_attributes=[('class', 'information')])
        self.report.add_html_object(section)

    def __add_overview_section(self) -> None:
        """
        Adds a new section which allows to download the generated assemblies at different stages of the pipeline.
        :return: None
        """
        section = HtmlReportSection('Overview section')
        rows = []
        for io_fasta in self._tool_inputs['FASTA']:
            fasta_key = io_fasta.path.parent.name
            relative_path = Path('assemblies', f'consensus_{fasta_key}.fasta')
            section.add_file(io_fasta.path, relative_path)
            rows.append([
                fasta_key, HtmlTableCell('Download (FASTA)', link=str(relative_path))
            ])

        # Add table with download links
        section.add_table(rows, column_names=['Assembly step', 'Download'], table_attributes=[('class', 'data')])
        section.add_paragraph(
            "The hybrid assembly pipeline performs long-read first <i>de novo</i> assembly, according to the following" 
            " steps: (1) Quality control and pre-processing of the long and short reads; (2) Long-reads assembly using "
            "Flye; (3) polishing using Medaka (long-reads), followed by Polypolish and Pypolca (short-reads); (4) "
            "Quality assessment using QUAST and several variant callers. An additional short-read first assembly is "
            "created using Unicycler (optional).")
        section.copy_files(self.report.output_dir)
        self.report.add_html_object(section)

    def __add_overview_links(self) -> None:
        """
        Adds the report header.
        :return: None
        """
        self.report.add_module_header('Sections')
        section = HtmlReportSection(None)
        overview_list = HtmlElement('ul')
        for title, key in HybridAssemblyReporter.REPORT_STRUCTURE:
            list_item = HtmlElement('li')
            list_item.add_html_object(HtmlElement('a', title, [('href', f'#{key}')]))
            overview_list.add_html_object(list_item)
        section.add_html_object(overview_list)
        self.report.add_html_object(section)

    def __add_quast_section(self, path_tsv: Path, quast_informs: dict[str, Any]) -> None:
        """
        Adds the summary QUAST table to the report.
        :param path_tsv: Path to the combined QUAST TSV file
        :param quast_informs: QUAST informs
        :return: None
        """
        section = HtmlReportSection('QUAST statistics', subtitle=quast_informs['_name_full'])
        data_quast = pd.read_table(path_tsv)
        section.add_table(
            list(data_quast.itertuples(index=False, name=None)),
            column_names=list(data_quast.columns),
            table_attributes=[('class', 'data')])
        relative_path = Path('quast', 'quast_all.html')
        section.add_file(self._tool_inputs['HTML_quast'][0].path, relative_path)
        section.add_link_to_file('Combined QUAST report (HTML)', relative_path)
        self.report.add_html_object(section)
        section.copy_files(self.report.output_dir)

    def __add_vc_table(self, path_tsv: Path) -> None:
        """
        Adds the variant calling section to the report.
        :path: Path to the tsv containing variant calling statistics
        :return: None
        """
        subtitle = ', '.join([self._input_informs[x]['_name_full'] for x in ('freebayes', 'clair3')])
        section = HtmlReportSection('Variant calling', subtitle=subtitle)
        data_vc = pd.read_table(path_tsv)
        section.add_table(
            list(data_vc.itertuples(index=False, name=None)),
            column_names=list(data_vc.columns),
            table_attributes=[('class', 'data')])
        section.add_paragraph(
            'Freebayes and Clair3 call SNPs and (small) indels on the short- and long-reads, respectively.')
        self.report.add_html_object(section)

    def __add_sniffles_table(self, path_tsv: Path) -> None:
        """
        Adds the Sniffles table to the report.
        :path: Path to the tsv containing variant calling from sniffles
        :return: None
        """
        section = HtmlReportSection('Sniffles statistics', subtitle=self._input_informs['sniffles'][0]['_name_full'])
        data_sniffles = pd.read_table(path_tsv)
        vcf_files = []
        for io_vcf in self._tool_inputs['VCF_sniffles']:
            vcf_key = io_vcf.path.parents[1].name
            relative_path = Path('sniffles', vcf_key, f'sniffles_{vcf_key}.vcf')
            vcf_files.append(HtmlTableCell('Download (VCF)', link=str(relative_path)))
            section.add_file(io_vcf.path, relative_path)
        data_sniffles['Download (VCF)'] = vcf_files
        section.add_table(
            list(data_sniffles.itertuples(index=False, name=None)),
            column_names=list(data_sniffles.columns),
            table_attributes=[('class', 'data')])
        section.add_paragraph("""
            Sniffles detect structural variation by mapping the long reads to the consensus sequence. Long deletions are 
            listed in the indels column, and all other variants (inversions, duplications, breakpoints) are listed in 
            the structural variation column (SV).""")
        self.report.add_html_object(section)
        section.copy_files(self.report.output_dir)

    def __add_mapping_table(self, path_tsv: Path) -> None:
        """
        Adds the mapping statistics to the report.
        :path_tsv: Path to the tsv containing mapping statistics
        :return: None
        """
        section = HtmlReportSection('Mapping statistics', subtitle=self._input_informs['mapping']['_name_full'])
        data_mapping = pd.read_table(path_tsv)
        section.add_table(
            list(data_mapping.itertuples(index=False, name=None)),
            column_names=list(data_mapping.columns),
            table_attributes=[('class', 'data')]
        )
        self.report.add_html_object(section)

    def __add_ale_score_table(self) -> None:
        """
        Adds the ALE score table, as well as downloadable Wiggle files to the report.
        :return: None
        """
        section = HtmlReportSection('ALE scores', subtitle=self._input_informs['ale'][0]['_name_full'])
        output_rows = []

        # Group input files
        wiggle_by_ale_key_by_assembly_key = {}
        for io_wiggle in self._tool_inputs['WIGGLE_ale']:
            ale_key = io_wiggle.path.stem.split('-')[-1]
            assembly_key = io_wiggle.path.parents[2].name
            if assembly_key not in wiggle_by_ale_key_by_assembly_key:
                wiggle_by_ale_key_by_assembly_key[assembly_key] = {}
            relative_path = Path('ale', f'{assembly_key}_{ale_key}.wig')
            wiggle_by_ale_key_by_assembly_key[assembly_key][ale_key] = relative_path
            section.add_file(io_wiggle.path, relative_path)

        # Collect scores
        score_by_assembly_key = {}
        for inform in self._input_informs['ale']:
            score_by_assembly_key[inform['_tag']] = inform['ale_score']

        # Create output table
        for assembly_key, wiggle_by_ale in wiggle_by_ale_key_by_assembly_key.items():
            output_rows.append([
                assembly_key,
                f'{int(score_by_assembly_key[assembly_key]):,}',
                *[HtmlTableCell('Download (WIG)', link=str(wiggle)) for ale_key, wiggle in wiggle_by_ale.items()]
            ])
        section.add_table(
            output_rows,
            column_names=[
                'Assembly step', 'Ale score', 'Download depth', 'Download K-mer', 'Download insert', 'Download place'],
            table_attributes=[('class', 'data')]
        )
        section.add_paragraph("""
            <b>Note:</b> The ALE score provides an indication of the overall quality of an assembly. Larger ALE scores 
            are better (ALE scores are always negative, which means that values closer to 0 are better)""")
        section.add_header('WIGGLE files', 4)
        section.add_table([
            ['depth', 'Agreement between the observed and expected depth, considering the %GC-content'],
            ['insert', 'Agreement between the observed and expected read pairing'],
            ['kmer', 'Likelihood of the assembly formula, in the absence of any read information'],
            ['place', 'Agreement between the sequences of reads and the assembly'],
        ], ['File', 'Explanation'], [('class', 'data')])
        section.add_paragraph(
            "WIGGLE files can be loaded into genome browser such as IGV to visualize issues with the assembly. More "
            "information is provided in the ALE manuscript, which is listed below (Clark et al., 2013).")
        self.report.add_html_object(section)
        section.copy_files(self.report.output_dir)
