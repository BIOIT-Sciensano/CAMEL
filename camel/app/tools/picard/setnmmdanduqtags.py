from camel.app.tools.picard.picard import Picard


class SetNmMdAndUqTags(Picard):
    """
    ==============================
    Picard SetNmMdAndUqTags 2.23.3
    ==============================
    Takes in a coordinate-sorted SAM or BAM and calculates the NM, MD, and UQ tags by comparing with the reference.

    Required inputs:
    ----------------
    ['BAM'|'SAM"]:      ToolIOFile object. A sorted SAM or BAM file.
    'FASTA_REF':        ToolIOFile object. FASTA file containing the reference genome.

    Output:
    -------
    'BAM'              ToolIOFile object. Fixed BAM file
    """

    def __init__(self):
        """
        Initialize a picard tool
        :return: None
        """
        super().__init__('Picard SetNmMdAndUqTags ', '2.23.3')

        self._required_inputs = ["BAM", "SAM", "FASTA_REF"]
