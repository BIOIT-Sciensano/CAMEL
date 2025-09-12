from camel.app.loggers import logger
from itertools import chain


class Pipeline:
    """
    Class meant to handle the workflow of steps.
    TODO: Deprecate?
    """

    def __init__(self, name: str, keep_config: bool = False, keep_error_log: bool = False) -> None:
        """
        Initializes a pipeline.
        :param name: Pipeline name
        :param keep_config: If True, pipeline config is saved on disk
        :param keep_error_log: If True, log is saved in case of error
        """
        self._name = name
        self._initial_input = None
        self._job_id = None
        self._keep_config = keep_config
        self._keep_error_log = keep_error_log
        logger.info(f"Initialized pipeline '{self._name}'")

    @property
    def name(self) -> str:
        """
        Returns the pipeline name.
        :return: Name
        """
        return self._name

    @property
    def keep_config(self) -> bool:
        """
        Returns True if the config file needs to be kept, False otherwise.
        :return: True if config is kept
        """
        return self._keep_config

    @property
    def keep_error_log(self) -> bool:
        """
        Returns True if the error log needs to be kept.
        :return: True if error log is kept
        """
        return self._keep_error_log

    def set_initial_input(self, files: dict) -> None:
        """
        Sets the initial inputs for the pipeline allowing it to be logged if logging is requested.
        :param files: Dictionary of input files
        :return: None
        """
        self._initial_input = files
        self._log_initial_input()

    def get_initial_input(self, key: str = None) -> list:
        """
        Returns a list of file paths so that they can be used in Snakemake. If a key is given only the file paths
        for that key are returned, otherwise all file paths in the dictionary are returned.
        :param key: Optional input file key
        :return: List of input file paths
        """
        io_files = list(chain(*self._initial_input.values())) if key is None else self._initial_input[key]
        return [file_.path for file_ in io_files]

    def _log_initial_input(self) -> None:
        """
        Logs the initial input of the pipeline.
        :return: None
        """
        for key, files in self._initial_input.items():
            for index, io in enumerate(files):
                if not io.is_logged:
                    continue
                logger.info(f"Pipeline input log: type={io.type_name} key={key} index={index} hash={io.hash}")
