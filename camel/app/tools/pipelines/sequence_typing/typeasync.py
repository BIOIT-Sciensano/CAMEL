import concurrent.futures
import json
from typing import Any, Union

from camel.app.toolkits.sequencetyping import typingasyncutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.errors import ToolExecutionError
from camel.app.loggers import logger
from camel.app.core.tool import Tool


class TypeAsync(Tool):
    """
    Performs BLAST-based sequence typing asynchronously for all loci of a scheme using a ThreadPool.
    This tool can be used to reduce the overhead in Snakemake by reducing the number of rules that need to be executed.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Typing async', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        if (self._parameters['detection_method'].value == 'blast') and ('FASTA' not in self._tool_inputs):
            raise InvalidToolInputError('FASTA input is required')
        elif (self._parameters['detection_method'].value in ('kma',)) and (not any(
                k in self._tool_inputs for k in ('FASTQ_SE', 'FASTQ_PE'))):
            raise InvalidToolInputError('FASTQ input is required')
        if 'DIR' not in self._tool_inputs:
            raise InvalidToolInputError("Typing directory input is required ('DIR')")
        super()._check_input()

    def __get_typing_parameters(self, locus_info: dict[str, Union[str, None]]) -> dict[str, Any]:
        """
        Returns the parameters to pass to the typing jobs.
        """
        # General options
        dict_params: dict[str, Any] = {
            'dir_working': self.folder / locus_info['name_sanitized'],
            'dir_scheme': self._tool_inputs['DIR'][0].path,
            'locus_metadata': locus_info,
            'threads_per_job': 1
        }
        # Add tool-specific options
        if self._parameters['detection_method'].value == 'blast':
            dict_params['fasta_in'] = self._tool_inputs['FASTA'][0].path
            dict_params['blastn_task'] = self._parameters['blastn_task'].value
        elif self._parameters['detection_method'].value == 'kma':
            dict_params['fastq_in'] = {k: v for k, v in self._tool_inputs.items() if k.startswith('FASTQ')}
        return dict_params

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        """
        # Parse scheme metadata
        with (self._tool_inputs['DIR'][0].path / 'scheme_metadata.txt').open() as handle:
            metadata = json.load(handle)
        logger.info(f"{len(metadata['loci'])} loci found")

        # Determine which loci need to be typed
        loci_target = [
            info['name_sanitized'] for info in metadata['loci'] if info['type'] == self._parameters['locus_type'].value]
        if len(loci_target) == 0:
            self._tool_outputs['VAL_hits'] = []
            return

        # Execute jobs in a thread pool
        output_by_locus = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=int(self._parameters['threads'].value)) as executor:
            detection_func = typingasyncutils.detection_by_method[self.get_param_value('detection_method')]
            future_to_locus = {
                executor.submit(detection_func, **self.__get_typing_parameters(locus_info)):
                    locus_info['name_sanitized'] for locus_info in metadata['loci'] if
                locus_info['name_sanitized'] in loci_target}

            for future in concurrent.futures.as_completed(future_to_locus):
                locus = future_to_locus[future]
                try:
                    data = future.result()
                except Exception as err:
                    logger.error(f"Locus '{locus}' generated exception: {err}")
                    raise ToolExecutionError(self.name, str(err))
                else:
                    output_by_locus[locus] = data

        # Set output
        self._tool_outputs['VAL_hits'] = [output_by_locus[locus].hit for locus in loci_target]

        # Collect informs (for the first locus)
        for key, value in output_by_locus[loci_target[0]].informs.items():
            self._informs[key] = value
