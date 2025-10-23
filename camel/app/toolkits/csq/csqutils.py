from dataclasses import dataclass


@dataclass
class BCSQInfo:
    """
    This class represents the info in the BCSQ info tag.
    """
    type_: str
    raw_str: str

    @staticmethod
    def parse(info_str: str) -> 'BCSQInfo':
        """
        Parses the BCSQ info tag.
        :param info_str: Info as a string
        :return: Parsed info
        """
        parts = info_str.split('|')
        return BCSQInfo(parts[0], info_str)
