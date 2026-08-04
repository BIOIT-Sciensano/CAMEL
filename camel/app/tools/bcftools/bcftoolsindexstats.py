from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class BcftoolsIndexStats(Tool):
    """
    Returns stats based on existing index files.
    Stats:
    - Number of variants per reference sequence ('variants_per_reference')
    - Total number of variants ('total_variants')
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('bcftools index stats', '1.17')

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        if not any(key in self._tool_inputs for key in ('BCF', 'VCF_GZ')):
            raise InvalidToolInputError("No input file found")
        if len(self._tool_inputs) != 1:
            raise InvalidToolInputError("Only one type of input is supported (VCF_GZ or BCF)")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        input_key = list(self._tool_inputs.keys())[0]
        self._command.command = f'{self._tool_command} {self._tool_inputs[input_key][0].path}'
        self._execute_command()
        self.__analyze_stdout()

    def __analyze_stdout(self) -> None:
        """
        Analyzes the standard output to report the SNP statistics.
        :return: None
        """
        self._informs['variants_per_reference'] = {}
        for line in self._command.stdout.splitlines():
            reference, size, snps = line.split('\t')
            self._informs['variants_per_reference'][reference] = int(snps)
        self._informs['total_variants'] = sum(
            self._informs['variants_per_reference'][ref] for ref in self._informs['variants_per_reference'].keys())
