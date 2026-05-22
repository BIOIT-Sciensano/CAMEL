from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlreportsection import HtmlReportSection

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class AMRCDSCompletenessReporter(Tool):
    """
    Class that checks the CDS completeness the genes associated with AMR.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('AMR CDS completeness', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'VAL_HITS' not in self._tool_inputs:
            raise InvalidToolInputError("'VAL_HITS' input is required")
        if 'DB' not in self._input_informs:
            raise InvalidToolInputError("'DB' informs input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the tool.
        :return: None
        """
        # Retrieve gene detection hits and group them by locus
        gene_detection_hits = [io.value for io in self._tool_inputs['VAL_HITS']]
        hit_by_locus = {h.locus: h for h in gene_detection_hits}

        # Get the minimum coverage parameter
        min_cov = float(self._parameters['min_cov'].value)
        self._informs['min_cov'] = min_cov

        # Determine which loci are missing
        mapping = self._input_informs['DB']['mapping']
        missing_loci = []
        self._informs['missing_loci'] = []
        for seq_id in mapping.keys():
            locus = mapping.get(seq_id)

            # Locus is present -> skip
            if (locus in hit_by_locus) and hit_by_locus[locus].subject_coverage >= min_cov:
                continue

            # Locus is present but below the threshold
            if (locus in hit_by_locus) and hit_by_locus[locus].subject_coverage < min_cov:
                missing_loci.append((
                    f'<i>{locus}</i>',
                    mapping.get_metadata(seq_id, 'antibiotics'),
                    f'{hit_by_locus[locus].subject_coverage:.2f}%'
                ))
                self._informs['missing_loci'].append(locus)
                continue

            # Locus is absent (i.e., below <20% coverage)
            missing_loci.append((f'<i>{locus}</i>', mapping.get_metadata(seq_id, 'antibiotics'), '<20%'))
            self._informs['missing_loci'].append(locus)

        # Create the output report section
        section = HtmlReportSection('CDS completeness')
        nb_targets = len(mapping)
        section.add_paragraph(
            f'{nb_targets - len(missing_loci)}/{nb_targets} coding regions in the resistance regions BED '
            f'file are covered at >{int(min_cov)}% with >90% identity (screening with KMA).')

        # Add table with missing loci
        if len(missing_loci) > 0:
            section.add_header('Missing loci', 4)
            section.add_table(
                [locus for locus in sorted(missing_loci)],
                ['Locus', 'Antibiotic(s)', 'Coverage'],
                [('class', 'data')])
            section.add_paragraph(
                'These CDSs are (partially) missing, and this must be taken into account when interpreting the AMR results.')

        # Set output
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]
