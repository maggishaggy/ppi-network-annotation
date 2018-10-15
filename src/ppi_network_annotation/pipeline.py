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
        gene_expression_file_path: str,
        entrez_id_header: str,
        log_fold_change_header: str,
        adjusted_p_value_header: str,
        split_char: str,
        maximum_adjusted_p_value: float,
        maximum_log2_fold_change: float,
        minimum_log2_fold_change: float,
        base_mean_header: Optional[str] = None,
        hippie_min_edge_weight: Optional[float] = None,
        current_disease_ids_path: Optional[str] = None,
        disease_associations_path: Optional[str] = None,
) -> Network:
    """Generate the protein-protein interaction network.

    :return Network: Protein-protein interaction network with information on differential expression.
    """
    # Compilation of a protein-protein interaction (PPI) graph (HIPPIE)
    protein_interactions = parsers.parse_ppi_graph(ppi_graph_path, hippie_min_edge_weight)
    protein_interactions = protein_interactions.simplify()

    if gene_expression_file_path.endswith('.xlsx'):
        gene_list = parsers.parse_excel(
            gene_expression_file_path,
            entrez_id_header=entrez_id_header,
            log_fold_change_header=log_fold_change_header,
            adjusted_p_value_header=adjusted_p_value_header,
            split_char=split_char,
            base_mean_header=base_mean_header,
        )
    elif gene_expression_file_path.endswith('.csv'):
        gene_list = parsers.parse_csv(
            gene_expression_file_path,
            entrez_id_header=entrez_id_header,
            log_fold_change_header=log_fold_change_header,
            adjusted_p_value_header=adjusted_p_value_header,
            split_char=split_char,
            base_mean_header=base_mean_header,
        )
    elif gene_expression_file_path.endswith('.tsv'):
        gene_list = parsers.parse_csv(
            gene_expression_file_path,
            entrez_id_header=entrez_id_header,
            log_fold_change_header=log_fold_change_header,
            adjusted_p_value_header=adjusted_p_value_header,
            split_char=split_char,
            base_mean_header=base_mean_header,
            sep="\t"
        )
    else:
        raise ValueError(f'Unsupported extension: {gene_expression_file_path}')

    if disease_associations_path is not None and current_disease_ids_path is not None:
        current_disease_ids = parsers.parse_disease_ids(current_disease_ids_path)
        disease_associations = parsers.parse_disease_associations(disease_associations_path,
                                                                  current_disease_ids)
    else:
        disease_associations = None

    # Build an undirected weighted graph with the remaining interactions based on Entrez gene IDs
    network = Network(
        protein_interactions,
        maximum_adjusted_p_value=maximum_adjusted_p_value,
        maximum_log2_fold_change=maximum_log2_fold_change,
        minimum_log2_fold_change=minimum_log2_fold_change,
    )
    network.set_up_network(gene_list, disease_associations=disease_associations)

    return network
