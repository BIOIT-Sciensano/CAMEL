import json
from pathlib import Path

from camel.app.components.mycobacterium import amrutils
from camel.app.components.mycobacterium.amrutils import ConfidenceLevel
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.tool import Tool
from camel.scripts.mycobacteriumpipeline import AMR_CIRCOS_TEMPLATE

#############
# Constants #
#############
MUTATION_COLORS = {
    ConfidenceLevel.ASSOC_R: 'vdred',
    ConfidenceLevel.ASSOC_R_int: 'vdred',
    ConfidenceLevel.ASSOC_S: 'vdgreen',
    ConfidenceLevel.ASSOC_S_int: 'vdgreen',
    ConfidenceLevel.UNKNOWN: 'vdgrey',
}

REGION_TYPE_COLORS = {
    'protein_coding': 'mbrown',
    'intergenic': 'morange',
    'pseudogene': 'mpink',
    'rRNA': 'myellow'
}

INNER_GENOME_R0 = 0.34
INNER_GENOME_WIDTH = 0.01

GENOME_COV_R0 = 0.22
GENOME_COV_WIDTH = 0.12

REGIONS_R0 = 0.48
REGIONS_WIDTH = 0.02

REGIONS_MUT_R0 = 0.73
REGIONS_MUT_WIDTH = 0.02


####################
# Class definition #
####################
class AMRCircosTemplateGeneration(Tool):
    """
    This tool creates the circos template to visualize the AMR mutations.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Mycobacterium: AMR circos template generation', '0.1')
        self._plots = []
        self._highlights = []
        self.__total_length = 0
        self._regions_by_name = {}
        self.__coord_by_region_name = {}
        self._mutations_by_idx = {}
        self.__coord_by_mut_idx = {}
        self.__config = []

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'TSV_depth' not in self._tool_inputs:
            raise InvalidToolInputError("TSV input with depth values is required ('TSV_depth')")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__calculate_coordinates(int(self._parameters['scale'].value), int(self._parameters['spacing'].value))
        karyogram = self.__create_karyogram()
        self.__add_inner_genome_circle(int(self._parameters['genome_size'].value))
        self.__add_genome_coverage()
        self.__add_regions()
        self.__add_regions_coverage()
        self.__add_connectors()
        self.__add_mutations()
        config_path = self.__export_config_file(karyogram)
        self._tool_outputs['TXT'] = [ToolIOFile(config_path)]

    def __calculate_coordinates(self, scale: int, spacing: int) -> None:
        """
        Calculates a novel coordinates system based on the regions that contain mutations.
        :param scale: Scale
        :param spacing: Spacing between regions
        :return: None
        """
        with open(self._tool_inputs['JSON'][0].path) as handle:
            mutation_data = json.load(handle)

        for idx, row in enumerate(mutation_data):
            if row['lofreq'] is True:
                continue
            associations = row['associations']
            if not any(amrutils.ConfidenceLevel(a['confidence']) in (
                    amrutils.ConfidenceLevel.ASSOC_R, amrutils.ConfidenceLevel.ASSOC_R_int) for a in associations):
                continue
            self._mutations_by_idx[idx] = row
            if row['region']['locus'] not in self._regions_by_name:
                self._regions_by_name[row['region']['locus']] = row['region']

        # Calculate region positions in the novel coordinate system
        current_pos = 0
        for region in sorted(self._regions_by_name.values(), key=lambda reg: reg['start']):
            region_size = region['end'] - region['start']
            new_start = current_pos * scale
            new_end = (current_pos + region_size) * scale
            current_pos += region_size + spacing
            self.__coord_by_region_name[region['locus']] = (new_start, new_end)
        self._total_size = current_pos * scale

        # If there are no regions, use a fixed total size
        if len(self._regions_by_name) == 0:
            self._total_size = 1000 * scale
        logger.info(f"Total size: {self._total_size}")

        # Calculate mutation positions in the novel coordinate system
        for idx, row in self._mutations_by_idx.items():
            new_start = self.__coord_by_region_name[row['region']['locus']][0]
            relative_pos = row['position'] - row['region']['start']
            new_position = new_start + relative_pos * scale
            self.__coord_by_mut_idx[idx] = new_position

    def __create_karyogram(self) -> Path:
        """
        Creates the karyogram.
        :return: Karyogram file path
        """
        karyogram_path = self._folder / 'karyogram.txt'
        logger.info("Creating karyogram data file: {}".format(karyogram_path))
        with open(karyogram_path, 'w') as handle_out:
            handle_out.write('chr - myco1 H37Rv 0 {} chr'.format(self._total_size))
            handle_out.write('\n')
        return karyogram_path

    def __add_inner_genome_circle(self, genome_size: int) -> None:
        """
        Adds a visualization of the inner genome.
        :param genome_size: genome size (in bp)
        :return: None
        """
        # Full circle
        path_full_circle = self._folder / 'highlights-full.txt'
        logger.info("Creating full circle: {}".format(path_full_circle))
        with open(path_full_circle, 'w') as handle:
            handle.write(' '.join(['myco1', '0', str(self._total_size)]))
            handle.write('\n')
        self._highlights.append(
            f"""
            <highlight>
            file = {path_full_circle}
            r1 = {INNER_GENOME_R0 + INNER_GENOME_WIDTH}r
            r0 = {INNER_GENOME_R0}r
            fill_color = black
            </highlight>""")
        self._plots.append(
            f"""
            <plot>
            type = highlight
            file = {path_full_circle}
            r1 = {GENOME_COV_R0 + 0.005}r
            r0 = {GENOME_COV_R0 - 0.005}r
            fill_color = black
            z = 1
            </highlight>""")

        # Small ticks
        interval = 1e5
        path_ticks_small = self._folder / 'ticks-small.txt'
        logger.info("Creating chromosome ticks: {}".format(path_ticks_small))
        with open(path_ticks_small, 'w') as handle:
            i = 0
            while i < genome_size:
                novel_pos = int((i / genome_size) * self._total_size)
                handle.write(' '.join(['myco1', str(novel_pos), str(novel_pos + 1)]))
                handle.write('\n')
                i += interval
        self._highlights.append(
            f"""
            <highlight>
            file = {path_ticks_small}
            r1 = {INNER_GENOME_R0 + INNER_GENOME_WIDTH + 0.005}r
            r0 = {INNER_GENOME_R0 + INNER_GENOME_WIDTH}r

            stroke_color=black
            stroke_thickness = 8
            z = 10
            </highlight>""")

        # Large ticks
        interval = 1e6
        tick_locations = []
        path_ticks_large = self._folder / 'ticks-large.txt'
        logger.info("Creating chromosome ticks: {}".format(path_ticks_large))
        with open(path_ticks_large, 'w') as handle:
            i = 0
            while i < genome_size:
                novel_pos = int((i / genome_size) * self._total_size)
                tick_locations.append(novel_pos)
                handle.write(' '.join(['myco1', str(novel_pos), str(novel_pos + 1)]))
                handle.write('\n')
                i += interval

        self._highlights.append(
            f"""
            <highlight>
            file = {path_ticks_large}
            r1 = {INNER_GENOME_R0 + INNER_GENOME_WIDTH + 0.02}r
            r0 = {INNER_GENOME_R0 + INNER_GENOME_WIDTH}r

            stroke_color=black
            stroke_thickness = 8
            z = 10
            </highlight>""")

        # Labels for the large ticks
        path_tick_labels = self.folder / 'labels-ticks.txt'
        logger.info("Creating tick labels: {}".format(path_tick_labels))
        with open(path_tick_labels, 'w') as handle:
            for i in range(0, len(tick_locations)):
                handle.write(' '.join(['myco1', str(tick_locations[i]), str(tick_locations[i] + 1), f'{i}Mb']))
                handle.write('\n')
        self._plots.append(
            f"""
            <plot>
            type = text
            file = {path_tick_labels}
            r1 = {INNER_GENOME_R0 + 0.20}r
            r0 = {INNER_GENOME_R0 + 0.04}r

            color = black
            z = 10
            label_size = 32
            </plot>""")

    def __add_genome_coverage(self) -> None:
        """
        Adds a plot with the genome coverage.
        :return: None
        """
        coverage_values = self.__parse_coverage_values(self._tool_inputs['TSV_depth'][0].path)
        step_size_scaled = self._total_size / max(len(coverage_values), 1)
        histogram_path = self.folder / 'hist-genome_cov.txt'
        with open(histogram_path, 'w') as handle:
            current_pos = 0
            for value in coverage_values:
                handle.write(
                    ' '.join(['myco1', str(int(current_pos)), str(int(current_pos + step_size_scaled)), str(value)])
                )
                handle.write('\n')
                current_pos += step_size_scaled

        self._plots.append(
            f"""
            <plot>
                type = histogram
                file = {histogram_path}
                r1 = {GENOME_COV_R0 + GENOME_COV_WIDTH}r
                r0 = {GENOME_COV_R0}r

                extend_bin = no
                color = sciensano_dgreen
                fill_color = sciensano_dgreen
                z = 4
                min = 0
                max = 100
                <axes>
                    <axis>
                        z = 5
                        color = grey
                        thickness = 1
                        spacing = 0.10r
                    </axis>
                </axes>
                <backgrounds>
                    <background>
                        y1 = 1r
                        y0 = 0r
                        color = sciensano_lgreen
                    </background>
                </backgrounds>
            </plot>
            """
        )

    def __parse_coverage_values(self, tsv_file: Path, window_size: int = 500) -> list[int]:
        """
        Parses the coverage values reported by samtools depth.
        :param tsv_file: TSV file containing depth values
        :param window_size: Window size
        :return: List of coverage values
        """
        coverage_values = []
        buffer = []
        with open(tsv_file, 'r') as handle:
            for line in handle.readlines():
                buffer.append(int(line.strip().split('\t')[-1]))
                if len(buffer) > window_size:
                    coverage_values.append(int(sum(buffer) / len(buffer)))
                    buffer.clear()
        logger.info(f"length cov values: {len(coverage_values)} (window size: {window_size})")
        return coverage_values

    def __export_config_file(self, karyogram_path: Path) -> Path:
        """
        Exports the config file.
        :param karyogram_path: Karyogram path
        :return: Path to config file
        """
        circos_config_path = self._folder / 'circos-config.txt'
        logger.info("Creating config file for circos: {}".format(circos_config_path))
        with open(self._folder / 'circos-config.txt', 'w') as handle_out:
            with open(AMR_CIRCOS_TEMPLATE) as handle_in:
                handle_out.write(handle_in.read().format(
                    karyo=karyogram_path,
                    highlights='\n'.join(self._highlights),
                    plots='\n'.join(self._plots)
                ))
        return circos_config_path

    def __add_regions(self) -> None:
        """
        Adds the visualization of the regions of interest.
        :return: None
        """
        # AMR associated regions, colored by type
        region_color_highlight_path = self.folder / 'highlight-regions-by-type.txt'
        logger.info("Creating gene highlights data file: {}".format(region_color_highlight_path))
        with region_color_highlight_path.open('w') as handle:
            for region_name, converted_coordinates in self.__coord_by_region_name.items():
                handle.write(' '.join([
                    'myco1', str(converted_coordinates[0]), str(converted_coordinates[1]),
                    'fill_color={}'.format(REGION_TYPE_COLORS[self._regions_by_name[region_name]['type']])
                ]))
                handle.write('\n')

        # Inner
        self._plots.append(
            f"""
            <plot>
            type = highlight
            file = {region_color_highlight_path}
            r1 = {REGIONS_R0 + REGIONS_WIDTH}r
            r0 = {REGIONS_R0}r
            stroke_color = black
            stroke_thickness = 4
            </plot>
            """)

        # Add boxes for regions
        region_highlight_path = self._folder / 'highlight-regions.txt'
        logger.info("Creating gene highlights data file: {}".format(region_highlight_path))
        with region_highlight_path.open('w') as handle:
            for region, converted_coordinates in self.__coord_by_region_name.items():
                handle.write(' '.join([
                    'myco1', str(converted_coordinates[0]), str(converted_coordinates[1])
                ]))
                handle.write('\n')

        # Fill
        self._highlights.append(
            f"""
            <highlight>
            file = {region_highlight_path}
            r1 = {REGIONS_MUT_R0 + REGIONS_MUT_WIDTH}r
            r0 = {REGIONS_R0}r
            fill_color = vvlblue
            stroke_thickness = 2
            stroke_color = dblue
            z = -2
            </highlight>
            """)

        # Add labels
        region_labels_path = self.folder / 'labels-region-names.txt'
        logger.info("Creating region name labels: {}".format(region_labels_path))
        with open(region_labels_path, 'w') as handle:
            for region_name, converted_coordinates in self.__coord_by_region_name.items():
                middle = int((converted_coordinates[0] + converted_coordinates[1]) / 2)
                offset = int(0.0085 * self._total_size)

                handle.write(' '.join([
                    'myco1',
                    str(middle),
                    str(middle - offset),
                    self._regions_by_name[region_name]['locus'],
                    'label_font=bolditalic,label_size=34'
                ]))
                handle.write('\n')
                handle.write(' '.join([
                    'myco1',
                    str(middle),
                    str(middle + offset),
                    f"[{self._regions_by_name[region_name]['abs_short'].replace(', ', ',')}]",
                    'label_font=bold,label_size=26'
                ]))
                handle.write('\n')

        self._plots.append(
            f"""
            <plot>
            type = text
            file = {region_labels_path}
            r1 = {REGIONS_R0 + REGIONS_WIDTH + 0.50}r
            r0 = {REGIONS_R0 + REGIONS_WIDTH}r
            color = black
            label_snuggle=yes
            rpadding = 18p
            </plot>
            """)

    def __add_mutations(self) -> None:
        """
        Adds the circle with the visualization of the mutations.
        :return: None
        """
        # Add boxes for regions
        region_highlight_path = self.folder / 'highlight-regions-muts.txt'
        logger.info("Creating gene highlights data file: {}".format(region_highlight_path))
        with open(region_highlight_path, 'w') as handle:
            for region, converted_coordinates in self.__coord_by_region_name.items():
                handle.write(' '.join([
                    'myco1', str(converted_coordinates[0]), str(converted_coordinates[1])
                ]))
                handle.write('\n')

        # Outer regions
        self._highlights.append(
            f"""
            <highlight>
            file = {region_highlight_path}
            r1 = {REGIONS_MUT_R0 + REGIONS_MUT_WIDTH}r
            r0 = {REGIONS_MUT_R0}r
            stroke_color = black
            stroke_thickness = 4
            fill_color=sciensano_lgreen
            </highlight>
            """)

        # Get the color for each mutation
        color_by_mut_idx = {}
        for mut_idx in self.__coord_by_mut_idx.keys():
            associations = self._mutations_by_idx[mut_idx]['associations']
            for level in ConfidenceLevel:
                if any(a['confidence'] == level.value for a in associations):
                    color_by_mut_idx[mut_idx] = MUTATION_COLORS[level]

        # Add the highlights on the outer bar
        mutation_highlights_path = self.folder / 'highlight-mutations.txt'
        logger.info(f"Creating mutation highlights: {mutation_highlights_path}")
        with open(mutation_highlights_path, 'w') as handle:
            for idx, converted_position in self.__coord_by_mut_idx.items():
                handle.write('\t'.join([
                    'myco1', str(converted_position), str(converted_position + 1),
                    'stroke_color={}'.format(color_by_mut_idx[idx])]))
                handle.write('\n')
        self._plots.append(
            f"""
            <plot>
            type = highlight
            file = {mutation_highlights_path}
            r1 = {REGIONS_MUT_R0 + REGIONS_MUT_WIDTH}r
            r0 = {REGIONS_MUT_R0}r
            stroke_thickness = 8
            z = 1
            </plot>
            """)

        # Add the connector lines
        self._highlights.append(
            f"""
            <highlight>
            file = {mutation_highlights_path}
            r1 = {REGIONS_MUT_R0}r
            r0 = {REGIONS_R0 + REGIONS_WIDTH}r
            fill_color = lgrey
            </plot>
            """
        )

        # Add labels
        mutation_labels_path = self.folder / 'labels-mutations.txt'
        logger.info("Creating mutation highlights: {}".format(mutation_labels_path))
        with open(mutation_labels_path, 'w') as handle:
            for mut_idx, converted_position in self.__coord_by_mut_idx.items():
                mutation = self._mutations_by_idx[mut_idx]
                handle.write('\t'.join([
                    'myco1', str(converted_position), str(converted_position + 1),
                    mutation['name'].replace(' ', '_'),
                    f'color={color_by_mut_idx[mut_idx]}']))
                handle.write('\n')
        self._plots.append(
            f"""
            <plot>
            type = text
            file = {mutation_labels_path}
            r1 = {REGIONS_MUT_R0 + 0.30}r
            r0 = {REGIONS_MUT_R0 + REGIONS_MUT_WIDTH}r
            label_size = 34
            label_snuggle=yes
            rpadding = 18p
            </plot>
            """
        )

    def __add_connectors(self) -> None:
        """
        Adds lines that connect the boxes to their location on the chromosome.
        """
        connectors_path = self.folder / 'connectors.txt'
        logger.info("Creating connectors: {}".format(connectors_path))
        with open(connectors_path, 'w') as handle:
            for r_name, coordinates in self.__coord_by_region_name.items():
                outer_mid_coord = int((coordinates[0] + coordinates[1]) / 2)
                original_mid_point = (
                    self._regions_by_name[r_name]['end'] + self._regions_by_name[r_name]['start']) / 2
                inner_mid_coord = int(self._total_size * (
                        original_mid_point / int(self._parameters['genome_size'].value)))
                handle.write(' '.join(['myco1', str(inner_mid_coord), str(outer_mid_coord)]))
                handle.write('\n')
        self._plots.append(
            f"""
            <plot>
            type = connector
            file = {connectors_path}
            r1 = {REGIONS_R0}r
            r0 = {INNER_GENOME_R0 + INNER_GENOME_WIDTH}r
            connector_dims=0,0,0.80,0.20,0
            thickness = 4
            color = grey
            </plot>
            """
        )

    def __add_regions_coverage(self) -> None:
        """
        Adds a plot that shows the coverage of the regions.
        :return: None
        """
        step_size = 10
        cov_values = self.__parse_coverage_values(self._tool_inputs['TSV_depth'][0].path, step_size)
        region_coverage_path = self._folder / 'histogram-region-coverage.txt'
        logger.info("Creating gene highlights data file: {}".format(region_coverage_path))
        with open(region_coverage_path, 'w') as handle:
            for region_name, converted_coordinates in self.__coord_by_region_name.items():
                region = self._regions_by_name[region_name]
                start_block_index = int(region['start'] / (step_size + 1))
                nb_blocks = int(abs(region['start'] - region['end']) / (step_size + 1))
                scaled_step_size = abs(converted_coordinates[1] - converted_coordinates[0]) / nb_blocks
                for i in range(0, nb_blocks):
                    try:
                        handle.write(' '.join([
                            'myco1',
                            str(int(converted_coordinates[0] + (i * scaled_step_size))),
                            str(int(converted_coordinates[0] + ((i + 1) * scaled_step_size))),
                            str(int(cov_values[start_block_index + i]))
                        ]))
                        handle.write('\n')
                    except IndexError:
                        logger.warning(f"Index error: {start_block_index + i}, cov values size: {len(cov_values)}")
        self._plots.append(
            f"""
            <plot>
            type = histogram
            file = {region_coverage_path}
            r1 = {REGIONS_MUT_R0 - 0.01}r
            r0 = {REGIONS_R0 + REGIONS_WIDTH}r
            extend_bin = no
            color = vlblue
            fill_color = lblue
            z = -1
            </plot>
            """)
