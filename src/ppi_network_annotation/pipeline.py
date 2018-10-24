#!/usr/bin/env python3

"""
"""

import logging
import os
from typing import Optional

from ppi_network_annotation import parsers
from ppi_network_annotation.model.network import Network

__all__ = [
    'generate_ppi_network'
]

logger = logging.getLogger(__name__)

HERE = os.path.abspath(os.path.dirname(__file__))


def generate_ppi_network(
        ppi_graph_path: str,
        dge_path: str,
        entrez_id_header: str,
        log2_fold_change_header: str,
        adj_p_header: str,
        entrez_delimiter: str,
        max_adj_p: float,
        max_log2_fold_change: float,
        min_log2_fold_change: float,
        base_mean_header: Optional[str] = None,
        ppi_edge_min_confidence: Optional[float] = None,
        current_disease_ids_path: Optional[str] = None,
        disease_associations_path: Optional[str] = None,
) -> Network:
    """Generate the protein-protein interaction network.

    :return Network: Protein-protein interaction network with information on differential expression.
    """
    # Compilation of a protein-protein interaction (PPI) graph (HIPPIE)
    protein_interactions = parsers.parse_ppi_graph(ppi_graph_path, ppi_edge_min_confidence)
    protein_interactions = protein_interactions.simplify()

    if dge_path.endswith('.xlsx'):
        gene_list = parsers.parse_excel(
            dge_path,
            entrez_id_header=entrez_id_header,
            log_fold_change_header=log2_fold_change_header,
            adjusted_p_value_header=adj_p_header,
            split_char=entrez_delimiter,
            base_mean_header=base_mean_header,
        )
    elif dge_path.endswith('.csv'):
        gene_list = parsers.parse_csv(
            dge_path,
            entrez_id_header=entrez_id_header,
            log_fold_change_header=log2_fold_change_header,
            adjusted_p_value_header=adj_p_header,
            split_char=entrez_delimiter,
            base_mean_header=base_mean_header,
        )
    elif dge_path.endswith('.tsv'):
        gene_list = parsers.parse_csv(
            dge_path,
            entrez_id_header=entrez_id_header,
            log_fold_change_header=log2_fold_change_header,
            adjusted_p_value_header=adj_p_header,
            split_char=entrez_delimiter,
            base_mean_header=base_mean_header,
            sep="\t"
        )
    else:
        raise ValueError(f'Unsupported extension: {dge_path}')

    if disease_associations_path is not None and current_disease_ids_path is not None:
        current_disease_ids = parsers.parse_disease_ids(current_disease_ids_path)
        disease_associations = parsers.parse_disease_associations(disease_associations_path,
                                                                  current_disease_ids)
    else:
        disease_associations = None

    # Build an undirected weighted graph with the remaining interactions based on Entrez gene IDs
    network = Network(
        protein_interactions,
        max_adj_p=max_adj_p,
        max_l2fc=max_log2_fold_change,
        min_l2fc=min_log2_fold_change,
    )
    network.set_up_network(gene_list, disease_associations=disease_associations)

    return network
