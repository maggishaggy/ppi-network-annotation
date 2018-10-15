# -*- coding: utf-8 -*-

"""This module contains the class Network."""

import logging
from collections import defaultdict
from typing import Dict, List, Optional

import numpy as np
from igraph import Graph

logger = logging.getLogger(__name__)


class Network:
    """Encapsulate proteins and their interactions in an igraph.Graph object."""

    def __init__(self, entire_ppi_graph: Graph, maximum_adjusted_p_value, maximum_log2_fold_change,
                 minimum_log2_fold_change):
        """Construct a Network object.

        :param Graph entire_ppi_graph: A network of protein interactions.
        :param config.Params cfp: An object that includes paths, cutoffs and
        other necessary information.
        """
        logger.info("Initializing Network")

        self.maximum_adjusted_p_value = maximum_adjusted_p_value
        self.maximum_log2_fold_change = maximum_log2_fold_change
        self.minimum_log2_fold_change = minimum_log2_fold_change
        # create deep copy of the graph graph
        self.graph = entire_ppi_graph.copy()

    def set_up_network(self, genes: List, gene_filter: bool = False,
                       disease_associations: Optional[Dict] = None):
        """Set up the network.

         Filter genes out if requested and add attributes to the vertices.

        :param genes: A list of Gene objects.
        :param gene_filter: Removes all genes that are not in list <genes> if True.
        :param disease_associations:
        """
        if gene_filter:
            self.filter_genes([gene.entrez_id for gene in genes])
        self._add_vertex_attributes(genes, disease_associations)
        self.print_summary("Graph of all genes")

    def filter_genes(self, relevant_entrez: list) -> None:
        """Filter out the genes that are not in list relevant_entrez.

        :param list relevant_entrez: Entrez IDs of genes which are to be kept.
        """
        logger.info("In filter_genes()")
        irrelevant_genes = self.graph.vs.select(name_notin=relevant_entrez)
        self.graph.delete_vertices(irrelevant_genes)

    def _add_vertex_attributes(self, genes: list, disease_associations: dict = None) -> None:
        """Add attributes to vertices.

        :param genes: A list of genes containing attribute information.
        """
        self._set_default_vertex_attributes()
        self._add_vertex_attributes_by_genes(genes)

        # compute up-regulated and down-regulated genes
        up_regulated = self.get_upregulated_genes()
        down_regulated = self.get_downregulated_genes()

        # set the attributes for up-regulated and down-regulated genes
        self.graph.vs(up_regulated.indices)["is_diff_expressed"] = True
        self.graph.vs(up_regulated.indices)["up_regulated"] = True
        self.graph.vs(down_regulated.indices)["is_diff_expressed"] = True
        self.graph.vs(down_regulated.indices)["down_regulated"] = True

        if disease_associations is not None:
            for target_id, disease_id_list in disease_associations.items():
                if target_id in self.graph.vs["name"]:
                    self.graph.vs.find(name=target_id)["associated_diseases"] = disease_id_list

        logger.info("Number of all differentially expressed genes is: {}".
                    format(len(up_regulated) + len(down_regulated)))

    def _set_default_vertex_attributes(self) -> None:
        """Assign default values on attributes to all vertices."""
        self.graph.vs["log2_fold_change"] = 0
        self.graph.vs["padj"] = 0.5
        self.graph.vs["symbol"] = self.graph.vs["name"]
        self.graph.vs["is_diff_expressed"] = False
        self.graph.vs["up_regulated"] = False
        self.graph.vs["down_regulated"] = False

    def _add_vertex_attributes_by_genes(self, genes: list) -> None:
        """Assign values to attributes on vertices.

        :param list genes: A list of Gene objects from which values will be extracted.
        """
        for gene in genes:
            try:
                vertex = self.graph.vs.find(name=str(gene.entrez_id)).index
                self.graph.vs[vertex]['log2_fold_change'] = gene.log2_fold_change
                self.graph.vs[vertex]['symbol'] = gene.symbol
                self.graph.vs[vertex]['padj'] = gene.padj
            except ValueError:
                pass

    def get_upregulated_genes(self):
        """Get genes that are up-regulated.

        :return: Up-regulated genes.
        """
        up_regulated = self.graph.vs.select(self._is_upregulated_gene)
        logger.info(
            "No. of up-regulated genes after laying on network: {}".format(len(up_regulated)))
        return up_regulated

    def _is_significantly_differentiated(self, v):
        return v.attributes()['padj'] < self.maximum_adjusted_p_value

    def _is_upregulated_gene(self, v):
        return (
            self._is_significantly_differentiated(v) and
            v.attributes()['log2_fold_change'] > self.minimum_log2_fold_change
        )

    def get_downregulated_genes(self):
        """Get genes that are down-regulated.

        :return: Down-regulated genes.
        """
        down_regulated = self.graph.vs.select(self._is_downregulated_gene)
        logger.info(
            "No. of down-regulated genes after laying on network: {}".format(len(down_regulated)))
        return down_regulated

    def _is_downregulated_gene(self, v):
        return (
            self._is_significantly_differentiated(v) and
            v.attributes()['log2_fold_change'] < self.maximum_log2_fold_change
        )

    def get_upregulated_genes_network(self) -> Graph:
        """Get the graph of up-regulated genes.

        :return Graph: Graph of up-regulated genes.
        """
        logger.info("In get_upregulated_genes_network()")

        deg_graph = self.graph.copy()  # deep copy graph
        not_diff_expr = self.graph.vs(up_regulated_eq=False)

        # delete genes which are not differentially expressed or have no connections to others
        deg_graph.delete_vertices(not_diff_expr.indices)
        deg_graph.delete_vertices(deg_graph.vs.select(_degree_eq=0))

        return deg_graph

    def get_downregulated_genes_network(self) -> Graph:
        """Get the graph of down-regulated genes.

        :return Graph: Graph of down-regulated genes.
        """
        logger.info("In get_downregulated_genes_network()")

        deg_graph = self.graph.copy()  # deep copy graph
        not_diff_expr = self.graph.vs(down_regulated_eq=False)

        # delete genes which are not differentially expressed or have no connections to others
        deg_graph.delete_vertices(not_diff_expr.indices)
        deg_graph.delete_vertices(deg_graph.vs.select(_degree_eq=0))

        return deg_graph

    def get_shortest_paths_graph(self, genes_to_keep: list = None,
                                 keep_isolated_nodes: bool = False):
        """Get the shortest paths graph between differentially expressed + special genes.

        :genes_to_keep list: A list of special genes.
        :keep_isolated_nodes bool: Removes the vertices with no neighbors when False.
        :return Graph: The shortest paths graph between special genes.
        """
        logger.info("In get_shortest_paths_graph()")
        sp_graph = self.graph.copy()
        weights = list(1 - np.array(sp_graph.es['weight']))
        sp_graph.es['weight'] = weights

        # Get the indices of all genes which are not to be deleted
        # (genes_to_keep + diff.expr.)
        relevant_gene_ind = self.graph.vs.select(is_diff_expressed=True).indices
        if genes_to_keep is not None:
            genes_to_keep_ind = self.graph.vs.select(name_in=genes_to_keep).indices
            relevant_gene_ind = set(relevant_gene_ind).union(set(genes_to_keep_ind))

        # Calculate the shortest paths between relevant genes and save the edges
        # that reside in shortest paths to a set
        shortest_path_edges = set()
        for ind in relevant_gene_ind:
            shortest_paths = sp_graph.get_shortest_paths(ind,
                                                         to=relevant_gene_ind,
                                                         weights='weight')
            for path in shortest_paths:
                for i in range(len(path) - 1):
                    eid = sp_graph.get_eid(path[i], path[i + 1])
                    shortest_path_edges.add(eid)

        # get and remove irrelevant edges
        irrelevant_edges = set(sp_graph.es.indices) - shortest_path_edges
        sp_graph.delete_edges(irrelevant_edges)
        if not keep_isolated_nodes:
            sp_graph.delete_vertices(sp_graph.vs.select(_degree_eq=0))

        return sp_graph

    def get_neighborhood_graph(self, node_name: str, order: int = 1) -> Graph:
        """Get the neighborhood graph of a node.

        :param str node_name: Node whose neighborhood graph is requested.
        :return Graph: Neighborhood graph
        """
        logger.info("In get_neighborhood_graph()")
        neighbors = list(self.get_neighbor_names(node_name, order))
        neighbor_network = self.graph.copy()
        neighbor_network.delete_vertices(
            self.graph.vs.select(name_notin=neighbors))
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

    def print_summary(self, heading: str) -> None:
        """Print the summary of a graph.

        :param str heading: Title of the graph.
        """
        logger.info(heading)
        logger.info("Number of nodes: {}".format(len(self.graph.vs)))
        logger.info("Number of edges: {}".format(len(self.graph.es)))

    def _get_differentially_expressed_genes(self, diff_type):
        """Get the differentially expressed genes based on diff_type.

        :param str diff_type: Differential expression type chosen by the user; all, down, or up.
        :return list: A list of differentially expressed genes.
        """
        if diff_type == "up":
            diff_expr = self.graph.vs.select(up_regulated_eq=True)
        elif diff_type == "down":
            diff_expr = self.graph.vs.select(down_regulated_eq=True)
        else:
            diff_expr = self.graph.vs.select(is_diff_expressed_eq=True)
        return diff_expr

    def get_neighborhood_overlap(self, node1, node2, connection_type=None):
        """ Get the intersection of two nodes's neighborhoods.

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

    def write_adj_list(self, path):
        adj_list = self.get_adjlist()

        with open(path, mode="w") as file:
            for i, line in enumerate(adj_list):
                print(i, *line, file=file)

    def get_adjlist(self) -> List[List[int]]:
        return self.graph.get_adjlist()

    def write_attribute_adj_list(self, path):
        """Write the bipartite attribute graph to a file.

        :param str path: Path to the output file.
        """
        att_mappings = self.get_attribute_mappings()

        with open(path, mode="w") as file:
            for k, v in att_mappings.items():
                print("{} {}".format(k, " ".join(str(e) for e in v)), file=file)

    def get_attribute_mappings(self):
        """Get a dictionary of mappings between vertices and enumerated attributes.

        :return: Dictionary of mappings between vertices and enumerated attributes.
        """
        att_ind_start = len(self.graph.vs)
        att_mappings = defaultdict(list)
        att_ind_end = self._add_differential_expression_attributes(att_ind_start, att_mappings)
        if "associated_diseases" in self.graph.vs.attributes():
            self._add_disease_association_attributes(att_ind_end, att_mappings)
        return att_mappings

    def _add_differential_expression_attributes(self, att_ind_start, att_mappings):
        """ Add differential expression information to the attribute mapping dictionary.

        :param int att_ind_start: Start index for enumerating the attributes.
        :param dict att_mappings: Dictionary of mappings between vertices and enumerated attributes.
        :return: End index for attribute enumeration.
        """
        up_regulated_ind = self.graph.vs.select(up_regulated_eq=True).indices
        down_regulated_ind = self.graph.vs.select(down_regulated_eq=True).indices
        rest_ind = self.graph.vs.select(is_diff_expressed_eq=False).indices

        self._add_attribute_values(att_ind_start + 1, att_mappings, up_regulated_ind)
        self._add_attribute_values(att_ind_start + 2, att_mappings, down_regulated_ind)
        self._add_attribute_values(att_ind_start + 3, att_mappings, rest_ind)
        return att_ind_start + 4

    def _add_attribute_values(self, value, att_mappings, indices):
        """Add an attribute value to the given vertices.

        :param int value: Attribute value.
        :param dict att_mappings: Dictionary of mappings between vertices and enumerated attributes.
        :param list indices: Indices of the vertices.
        """
        for i in indices:
            att_mappings[i].append(value)

    def _add_disease_association_attributes(self, att_ind_start, att_mappings):
        """Add disease association information to the attribute mapping dictionary.

        :param int att_ind_start: Start index for enumerating the attributes.
        :param dict att_mappings: Dictionary of mappings between vertices and enumerated attributes.
        """
        disease_mappings = self._get_disease_mappings(att_ind_start)
        for vertex in self.graph.vs:
            assoc_diseases = vertex["associated_diseases"]
            if assoc_diseases is not None:
                assoc_disease_ids = [disease_mappings[disease] for disease in assoc_diseases]
                att_mappings[vertex.index].extend(assoc_disease_ids)

    def _get_disease_mappings(self, att_ind_start):
        """Get a dictionary of enumerations for diseases.

        :param int att_ind_start: Starting index for enumeration.
        :return: Dictionary of disease, number pairs.
        """
        all_disease_ids = self._get_all_unique_diseases()
        disease_enum = enumerate(all_disease_ids, start=att_ind_start)
        disease_mappings = {}
        for num, dis in disease_enum:
            disease_mappings[dis] = num
        return disease_mappings

    def _get_all_unique_diseases(self):
        """Get all unique diseases that are known to the network.

        :return: All unique disease identifiers.
        """
        all_disease_ids = self.graph.vs["associated_diseases"]
        # remove None values from list
        all_disease_ids = [lst for lst in all_disease_ids if lst is not None]
        # flatten list of lists, get unique elements
        all_disease_ids = list(set([id for sublist in all_disease_ids for id in sublist]))
        return all_disease_ids

    def write_index_labels(self, drug_targets, output_path):
        """ Write the mappings between vertex indices and labels(drug target vs. not) to a file.

        :param list drug_targets: List of known drug targets.
        :param str output_path: Path to the output file.
        """
        label_mappings = self.get_index_labels(drug_targets)

        with open(output_path, "w") as file:
            for k, v in label_mappings.items():
                print(k, v, sep='\t', file=file)

    def get_index_labels(self, drug_targets):
        drug_target_ind = self.graph.vs.select(name_in=drug_targets).indices
        rest_ind = self.graph.vs.select(name_notin=drug_targets).indices
        label_mappings = {i: 1 for i in drug_target_ind}
        label_mappings.update({i: 0 for i in rest_ind})
        return label_mappings

    def get_attribute_from_indices(self, indices: list, attribute_name: str):
        return list(np.array(self.graph.vs[attribute_name])[indices])

    def to_json(self, graph, known_targets, pred_targets):
        # {
        #     "nodes": [
        #         {"id": "Myriel", "group": 1},
        #         {"id": "Mme.Hucheloup", "group": 8}
        #     ],
        #     "links": [
        #         {"source": "Napoleon", "target": "Myriel", "value": 1},
        #         {"source": "Mme.Hucheloup", "target": "Enjolras", "value": 1}
        #     ]
        # }
        pass
