from pathlib import Path

from camel.app.toolkits.blast.blasthitstatistics import BlastHitStatistics
from camel.app.toolkits.blasttyping import blasthitclustering
from camel.app.toolkits.blasttyping.blasthitfilteringhelper import BlastHitFilteringHelper
from camel.app.toolkits.genedetection.genedetectionblasthit import GeneDetectionBlastHit
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.errors import ToolExecutionError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.core.tool import Tool
from camel.app.toolkits.genedetection.genedetectionkmahit import GeneDetectionKMAHit

class MixedHitFiltering(Tool):
    """
    Performs mixed hit filtering (KMA + BLAST).
    """
    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Gene Detection: Mixed Hit Filtering', '0.1')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        pass

    def _check_input(self) -> None:
        """
        Checks whether the input is correct.
        :return: None
        """
        if 'TSV_KMA' not in self._tool_inputs:
            raise InvalidToolInputError("No 'TSV' input for KMA found.")
        if 'TSV_BLAST' not in self._tool_inputs:
            raise InvalidToolInputError("No 'TSV' input for Blast found.")
        super()._check_input()

    def __set_informs(self) -> None:
        """
        Sets the informs.
        :return: None
        """
        self._informs['Min. percent identity'] = f"{self._parameters['min_percent_identity'].value}%"
        self._informs['Min. coverage'] = f"{self._parameters['min_coverage'].value}%"

