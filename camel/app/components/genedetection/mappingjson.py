from typing import Any, Optional


class MappingJSON:
    """
    This class contains metadata for the standardized FASTA files of the gene detection workflow.
    The main purpose of this class is to avoid cluttering of the log with the complete mapping as a dictionary.
    """

    def __init__(self, content: dict[str, Any]) -> None:
        """
        Initializes the mapping.
        :param content: Mapping content
        :return: None
        """
        self._content = content

    def __repr__(self) -> str:
        """
        Returns the printable representation of the mapping.
        :return: Representation
        """
        return f'Mapping({len(self._content):,} items)'

    def __len__(self) -> int:
        """
        Returns the length of the mapping.
        """
        return len(self._content)

    def keys(self) -> list[str]:
        """
        Returns the keys of the mapping.
        """
        return list(self._content.keys())

    def get(self, seq_id: str) -> None:
        """
        Returns the original header.
        :param seq_id: Sequence id
        :return: Original header
        """
        return self._content[seq_id]['header_orig']

    def get_metadata(self, seq_id: str, metadata_key: str, default: Optional[str] = None) -> str:
        """
        Returns the metadata value for the given key and sequence identifier.
        :param seq_id: Sequence identifier
        :param metadata_key: Metadata key
        :param default: Default value if key not present in metadata
        :return: None
        """
        if seq_id not in self._content:
            raise ValueError(f"No sample with id '{seq_id}' in mapping")

        metadata = self._content[seq_id]
        if metadata_key not in metadata:
            if default is None:
                raise ValueError(f"Key '{metadata_key}' not found in metadata")
            else:
                return default
        return metadata[metadata_key]
