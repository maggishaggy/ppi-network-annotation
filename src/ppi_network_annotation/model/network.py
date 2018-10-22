# -*- coding: utf-8 -*-

"""This module contains the class Network."""

import logging
from typing import Dict, List, Optional

import numpy as np
from igraph import Graph

logger = logging.getLogger(__name__)


class Network:
    """Encapsulate a PPI network with differential gene expression and disease association annotation."""

    def __init__(self, ppi_graph: Graph, max_adj_p, max_l2fc, min_l2fc):
        """Initialize the network object

        :param ppi_graph: A graph of protein interactions.
        :param max_adj_p: Maximum value for adjusted p-value, used for calculating differential expression
        :param max_l2fc: Maximum value for log2 fold change, used for calculating down regulation
        :param min_l2fc: Minimum value for log2 fold change, used for calculating up regulation
        """
        logger.info("Initializing Network")

        self.max_adj_p = max_adj_p
        self.max_l2fc = max_l2fc
        self.min_l2fc = min_l2fc
        self.graph = ppi_graph.copy()  # create deep copy of the graph graph

    def set_up_network(self, genes: List, gene_filter: bool = False,
                       disease_associations: Optional[Dict] = None):
        """Set up the network.

         Filter genes out if requested and add attributes to the vertices.

        :param genes: A list of Gene objects.
        :param gene_filter: Removes all genes that are not in list <genes> if True.
        :param disease_associations: Diseases associated with genes.
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
        self.graph.vs(up_regulated.indices)["diff_expressed"] = True
        self.graph.vs(up_regulated.indices)["up_regulated"] = True
        self.graph.vs(down_regulated.indices)["diff_expressed"] = True
        self.graph.vs(down_regulated.indices)["down_regulated"] = True

        # add disease associations
        self._add_disease_associations(disease_associations)

        logger.info("Number of all differentially expressed genes is: {}".
                    format(len(up_regulated) + len(down_regulated)))

    def _set_default_vertex_attributes(self) -> None:
        """Assign default values on attributes to all vertices."""
        self.graph.vs["l2fc"] = 0
        self.graph.vs["padj"] = 0.5
        self.graph.vs["symbol"] = self.graph.vs["name"]
        self.graph.vs["diff_expressed"] = False
        self.graph.vs["up_regulated"] = False
        self.graph.vs["down_regulated"] = False

    def _add_vertex_attributes_by_genes(self, genes: list) -> None:
        """Assign values to attributes on vertices.

        :param list genes: A list of Gene objects from which values will be extracted.
        """
        for gene in genes:
            try:
                vertex = self.graph.vs.find(name=str(gene.entrez_id)).index
                self.graph.vs[vertex]['l2fc'] = gene.log2_fold_change
                self.graph.vs[vertex]['symbol'] = gene.symbol
                self.graph.vs[vertex]['padj'] = gene.padj
            except ValueError:
                pass

    def _add_disease_associations(self, disease_associations: dict) -> None:
        """Add disease association annotation to the network.

        :param disease_associations: Dictionary of disease-gene associations.
        """
        if disease_associations is not None:
            for target_id, disease_id_list in disease_associations.items():
                if target_id in self.graph.vs["name"]:
                    self.graph.vs.find(name=target_id)["associated_diseases"] = disease_id_list

    def get_upregulated_genes(self):
        """Get genes that are up-regulated.

        :return: Up-regulated genes.
        """
        up_regulated = self.graph.vs.select(self._is_upregulated_gene)
        logger.info(f"No. of up-regulated genes after laying on network: {len(up_regulated)}")
        return up_regulated

    def get_downregulated_genes(self):
        """Get genes that are down-regulated.

        :return: Down-regulated genes.
        """
        down_regulated = self.graph.vs.select(self._is_downregulated_gene)
        logger.info(f"No. of down-regulated genes after laying on network: {len(down_regulated)}")
        return down_regulated

    def _is_significantly_differentiated(self, v):
        return v.attributes()['padj'] < self.max_adj_p

    def _is_upregulated_gene(self, v):
        return self._is_significantly_differentiated(v) and v.attributes()['l2fc'] > self.min_l2fc

    def _is_downregulated_gene(self, v):
        return self._is_significantly_differentiated(v) and v.attributes()['l2fc'] < self.max_l2fc

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
            diff_expr = self.graph.vs.select(diff_expressed_eq=True)
        return diff_expr

    def write_adj_list(self, path):
        adj_list = self.get_adjlist()

        with open(path, mode="w") as file:
            for i, line in enumerate(adj_list):
                print(i, *line, file=file)

    def get_adjlist(self) -> List[List[int]]:
        return self.graph.get_adjlist()

    def get_disease_mappings(self, att_ind_start):
        """Get a dictionary of enumerations for diseases.

        :param int att_ind_start: Starting index for enumeration.
        :return: Dictionary of disease, number pairs.
        """
        all_disease_ids = self.get_all_unique_diseases()
        disease_enum = enumerate(all_disease_ids, start=att_ind_start)
        disease_mappings = {}
        for num, dis in disease_enum:
            disease_mappings[dis] = num
        return disease_mappings

    def get_all_unique_diseases(self):
        """Get all unique diseases that are known to the network.

        :return: All unique disease identifiers.
        """
        all_disease_ids = self.graph.vs["associated_diseases"]
        # remove None values from list
        all_disease_ids = [lst for lst in all_disease_ids if lst is not None]
        # flatten list of lists, get unique elements
        all_disease_ids = list(set([id for sublist in all_disease_ids for id in sublist]))
        return all_disease_ids

    def write_index_labels(self, targets, output_path):
        """ Write the mappings between vertex indices and labels(target vs. not) to a file.

        :param list targets: List of known targets.
        :param str output_path: Path to the output file.
        """
        label_mappings = self.get_index_labels(targets)

        with open(output_path, "w") as file:
            for k, v in label_mappings.items():
                print(k, v, sep='\t', file=file)

    def get_index_labels(self, targets):
        target_ind = self.graph.vs.select(name_in=targets).indices
        rest_ind = self.graph.vs.select(name_notin=targets).indices
        label_mappings = {i: 1 for i in target_ind}
        label_mappings.update({i: 0 for i in rest_ind})
        return label_mappings

    def get_attribute_from_indices(self, indices: list, attribute_name: str):
        return list(np.array(self.graph.vs[attribute_name])[indices])