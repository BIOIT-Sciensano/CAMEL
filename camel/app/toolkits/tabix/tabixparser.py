from pathlib import Path

from camel.app.loggers import logger

TABIX_ANNOTATION_INDICES = {
    'position': 1,
    'type': 18,
    'locus': 16,
    'mutation_pos': 25
}


def parse_tabix_annotation(annotation_file: Path) -> dict[int, dict[str, str]]:
    """
    Parses the file with the tabular annotation.
    :param annotation_file: Annotation file
    :return: Annotations (position: annotation)
    """
    logger.info(f"Parsing annotation file: {annotation_file}")
    annotations = {}
    with open(annotation_file) as handle:
        for line in handle.readlines():
            parts = line.strip().split('\t')
            position = int(parts[TABIX_ANNOTATION_INDICES['position']])
            annotations[position] = {k: parts[index] for k, index in TABIX_ANNOTATION_INDICES.items()}
    return annotations
