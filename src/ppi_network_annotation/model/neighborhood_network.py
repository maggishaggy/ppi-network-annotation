# -*- coding: utf-8 -*-

"""This module contains the class NeighborhoodNetwork."""

import logging

from igraph import Graph

from ppi_network_annotation.model.network import Network

logger = logging.getLogger(__name__)


class NeighborhoodNetwork:
    """Mimic encapsulation of a neighborhood network."""

    def __init__(self, network: Network):
        """Initialize the network object.

        :param network: A PPI network
        """
        logger.info("Initializing NeighborhoodNetwork")

        self.graph = network.graph

    def get_neighborhood_network(self, node_name: str, order: int = 1) -> Graph:
        """Get the neighborhood graph of a node.

        :param str node_name: Node whose neighborhood graph is requested.
        :return Graph: Neighborhood graph
        """
        logger.info("In get_neighborhood_graph()")
        neighbors = list(self.get_neighbor_names(node_name, order))
        neighbor_network = self.graph.copy()
        neighbor_network.delete_vertices(self.graph.vs.select(name_notin=neighbors))
        return neighbor_network

    def get_neighbor_names(self, node_name: str, order: int = 1) -> list:
        """Get the names of all neighbors of a node, and the node itself.

        :param node_name: Node whose neighbor names are requested.
        :return: A list of names of all neighbors of a node, and the node itself.
        """
        logger.info("In get_neighbor_names()")
        node = self.graph.vs.find(name=node_name)
        neighbors = self.graph.neighborhood(node, order=order)
        names = self.graph.vs[neighbors]["name"]
        names.append(node_name)
        return list(names)

    def get_neighborhood_overlap(self, node1, node2, connection_type=None):
        """Get the intersection of two nodes's neighborhoods.

        Neighborhood is defined by parameter connection_type.
        :param Vertex node1: First node.
        :param Vertex node2: Second node.
        :param Optional[str] connection_type: One of direct or second-degree. Defaults to direct.
        :return: Overlap of the nodes' neighborhoods.
        """
        if connection_type is None or connection_type == "direct":
            order = 1
        elif connection_type == "second-degree":
            order = 2
        else:
            raise Exception(
                "Invalid option: {}. Valid options are direct and second-degree".format(
                    connection_type)
            )

        neighbors1 = self.graph.neighborhood(node1, order=order)
        neighbors2 = self.graph.neighborhood(node2, order=order)
        return set(neighbors1).intersection(neighbors2)
