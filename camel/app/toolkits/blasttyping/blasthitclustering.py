from camel.app.toolkits.blasttyping.cluster import Cluster
from camel.app.toolkits.genedetection.genedetectionblasthit import GeneDetectionBlastHit
from camel.app.loggers import logger


class BlastHitClustering:
    """
    Class that clusters overlapping BlastHit objects.
    Hits are considered overlapping if:
    - They are located on the same contig
    - They have at least one base overlap
    """

    @staticmethod
    def cluster_overlapping(hits: list[GeneDetectionBlastHit]) -> list[Cluster]:
        """
        Clusters the given blast hits into groups of overlapping hits.
        :param hits: Hits
        :return: List of clusters
        """
        clusters = []
        while len(hits) > 0:
            new_cluster = Cluster(hits[0])
            for hit in hits[1:]:
                if new_cluster.overlaps(hit):
                    new_cluster.add_hit(hit)
            for hit in new_cluster.hits:
                hits.remove(hit)
            clusters.append(new_cluster)
        BlastHitClustering.__log_clusters(clusters)
        logger.info(f'{len(clusters)} clusters of overlapping hits found')
        return clusters

    @staticmethod
    def __log_clusters(clusters):
        """
        Logs the clusters.
        :param clusters: Clusters
        :return: None
        """
        for cluster in clusters:
            logger.debug(f'Cluster {clusters.index(cluster) + 1}')
            for hit in cluster.hits:
                logger.debug(str(hit))
