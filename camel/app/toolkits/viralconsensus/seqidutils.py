import re

PATTERN_SEQ_ID = re.compile(r'^[\w.]+-\w+')


def get_segment(seq_id: str) -> str:
    """
    Returns the segment of a sequence id (e.g. 'OQ839344.1-HA' -> 'HA').
    :param seq_id: Sequence id
    :return: Segment name
    """
    return seq_id.rsplit('-', 1)[-1]


def get_accession(seq_id: str) -> str:
    """
    Returns the identifier (accession) of a sequence id (e.g. 'OQ839344.1-HA' -> 'OQ839344.1'). Sequence ids without
    segment are returned unchanged.
    :param seq_id: Sequence id
    :return: Accession
    """
    return seq_id.rsplit('-', 1)[0]


def is_valid_seq_id(seq_id: str) -> bool:
    """
    Checks if the sequence id has the '{identifier}-{segment}' format.
    :param seq_id: Sequence id
    :return: True if valid
    """
    return PATTERN_SEQ_ID.match(seq_id) is not None
