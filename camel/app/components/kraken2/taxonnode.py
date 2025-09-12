import dataclasses
import typing


@dataclasses.dataclass(frozen=True, unsafe_hash=True)
class TaxonNode:
    """
    Helper class to parse the Kraken 2 report output as a taxonomic tree.
    """
    name: typing.Union[str, int]
    percentage: typing.Optional[float] = None
    children: set['TaxonNode'] = dataclasses.field(default_factory=set, compare=False, hash=False)

    def add_child(self, node: 'TaxonNode') -> None:
        """
        Adds a child node to this taxon node.
        :param node: node to add
        :return: None
        """
        if node not in self.children:
            self.children.add(node)

    def __repr__(self, level: int=0) -> None:
        """
        Returns the string representation of this taxon node.
        :param level: level of the node
        :return: string representation
        """
        indent = ' ' * (level * 2)
        repr_str = f"{indent}{self.name}\n"
        for child in self.children:
            repr_str += child.__repr__(level + 1)
        return repr_str

    def total_perc(self) -> float:
        """
        Calculates the total percentage of the taxon node.
        :return: Total percentage
        """
        base_percentage = self.percentage if self.percentage is not None else 0
        return base_percentage + sum(c.total_perc() for c in self.children)

    @staticmethod
    def traverse_dfs(node: 'TaxonNode') -> typing.Iterator['TaxonNode']:
        """
        Traverses this node and all its child nodes.
        :param node: node to traverse
        :return: Node
        """
        yield node
        for child in sorted(node.children, key=lambda x: x.name):
            yield from TaxonNode.traverse_dfs(child)
