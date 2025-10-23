import json

from camel.app.toolkits.sequencetyping.sequencetypingutils import LocusMetadataHolder
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.errors import ToolExecutionError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool


class TypingDBLoader(Tool):
    """
    This tool extracts the metadata from a scheme / locus set directory.

    Input:
        DIR: Directory with several subdirectories for each locus and a scheme_metadata.txt file.

    Informs:
        The scheme metadata is added directly to the informs.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Typing: Locus Set Manager', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the specified input is correct.
        :return: None
        """
        if 'DIR' not in self._tool_inputs:
            raise InvalidToolInputError("DIR input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        path_dir_scheme = self._tool_inputs['DIR'][0].path

        # Parse metadata
        path_metadata = path_dir_scheme / 'scheme_metadata.txt'
        if not path_metadata.exists():
            raise ToolExecutionError(self.name, f"No scheme metadata found in '{path_metadata}'")
        with open(path_metadata) as handle:
            metadata = json.load(handle)
            metadata['loci'] = LocusMetadataHolder(metadata['loci'])
            self._informs.update(metadata)

        # Parse updating information
        path_json_update = path_dir_scheme / 'db_update_info.json'
        if not path_json_update.exists():
            raise FileNotFoundError(f'No database update JSON file found: {path_json_update}')
        with path_json_update.open() as handle:
            self._informs['update'] = json.load(handle)

        # Add the profiles (if they exist)
        path_profiles = path_dir_scheme / 'profiles.tsv'
        if path_profiles.exists():
            self._tool_outputs['TSV'] = [ToolIOFile(path_profiles)]
