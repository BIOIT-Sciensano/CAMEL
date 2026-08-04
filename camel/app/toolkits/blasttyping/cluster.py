from camel.app.toolkits.genedetection.genedetectionblasthit import GeneDetectionBlastHit


class Cluster:
    """
    Class that represents a cluster of BlastHit objects.
    """

    def __init__(self, hit: GeneDetectionBlastHit) -> None:
        """
        Initializes a cluster.
        :param hit: Initial hit
        """
        self.hits = [hit]
        self.seq_id = hit.blast_stats.query_id
        self._region = set(range(hit.blast_stats.query_start, hit.blast_stats.query_end))

    def add_hit(self, hit: GeneDetectionBlastHit) -> None:
        """
        Adds a hit to the cluster.
        :param hit: Blast hit
        :return: None
        """
        self.hits.append(hit)
        self._region = self._region.union(list(range(hit.blast_stats.query_start, hit.blast_stats.query_end)))

    def overlaps(self, hit):
        """
        Checks if the given hit overlaps with this cluster.
        :param hit: Blast hit
        :return: True if the hit overlaps with this cluster
        """
        if hit.blast_stats.query_id != self.seq_id:
            return False
        # noinspection PyTypeChecker
        return len(self._region.intersection(list(range(hit.blast_stats.query_start, hit.blast_stats.query_end)))) != 0
