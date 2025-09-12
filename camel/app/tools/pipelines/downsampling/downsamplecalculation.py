import json
import statistics
from typing import Any

from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class DownsampleCalculation(Tool):
    """
    Calculates the factor for downsampling based on the estimated coverage.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Downsample calculation', '0.1')

    def _check_input(self) -> None:
        """
        Check if the provided tool input is valid.
        :return: None
        """
        if 'stats' not in self._input_informs:
            raise InvalidToolInputError("Stats input is required")
        super()._check_input()

    def __calculate_stats(self, fq_stats: list[dict[str, Any]]) -> dict[str, Any]:
        """
        Calculates the downsampling statistics.
        :param fq_stats: FASTQ stats.
        :return: statistics as a dictionary
        """
        key_nb_reads = 'nb_read_pairs' if 'is_paired' in self._parameters else 'nb_reads'
        data_out = {
            'total_bases': sum([fq['nb_of_bases'] for fq in fq_stats]),
            'mean_read_length': statistics.mean([fq['nb_of_bases'] / fq['nb_of_sequences'] for fq in fq_stats]),
            f'{key_nb_reads}_in': next(iter(fq_stats))['nb_of_sequences'],
        }
        # Reference genome size is unknown -> cov. calculation is not possible
        if self._parameters['size_ref_genome'].value is None:
            data_out['downsample_factor'] = None
            return data_out

        # Calculate coverage
        ref_genome_size = int(self._parameters['size_ref_genome'].value)
        cov_est = sum(fq['nb_of_bases'] for fq in fq_stats) / ref_genome_size
        cov_target = float(self._parameters['cov_target'].value)
        downsample_factor = float(f"{cov_target / cov_est:.6f}")
        data_out.update({
            'coverage_estimated': cov_est,
            'coverage_target': cov_target,
            'downsample_factor': cov_target / cov_est if downsample_factor < 1 else None,
            'size_ref_genome': ref_genome_size,
        })
        return data_out

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        data_stats = self.__calculate_stats(self._input_informs['stats']['stats'])
        self._informs['stats'] = data_stats

        # Create JSON output file
        path_json_out = self.folder / 'downsample_stats.json'
        with path_json_out.open('w') as handle:
            json.dump(data_stats, handle, indent=2)
        self._tool_outputs['JSON'] = [ToolIOFile(path_json_out)]
