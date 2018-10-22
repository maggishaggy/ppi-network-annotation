from collections import defaultdict
import logging
from ppi_network_annotation.model.network import Network

logger = logging.getLogger(__name__)


class AttributeNetwork:
    def __init__(self, network: Network):
        """Initialize the network object

        :param network: A PPI network annotated with differential gene expression and disease association.
        """
        logger.info("Initializing AttributeNetwork")

        self.network = network

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
        att_ind_start = len(self.network.graph.vs)
        att_mappings = defaultdict(list)
        att_ind_end = self._add_differential_expression_attributes(att_ind_start, att_mappings)
        if "associated_diseases" in self.network.graph.vs.attributes():
            self._add_disease_association_attributes(att_ind_end, att_mappings)
        return att_mappings

    def _add_differential_expression_attributes(self, att_ind_start, att_mappings):
        """ Add differential expression information to the attribute mapping dictionary.

        :param int att_ind_start: Start index for enumerating the attributes.
        :param dict att_mappings: Dictionary of mappings between vertices and enumerated attributes.
        :return: End index for attribute enumeration.
        """
        up_regulated_ind = self.network.graph.vs.select(up_regulated_eq=True).indices
        down_regulated_ind = self.network.graph.vs.select(down_regulated_eq=True).indices
        rest_ind = self.network.graph.vs.select(diff_expressed_eq=False).indices

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
        disease_mappings = self.network.get_disease_mappings(att_ind_start)
        for vertex in self.network.graph.vs:
            assoc_diseases = vertex["associated_diseases"]
            if assoc_diseases is not None:
                assoc_disease_ids = [disease_mappings[disease] for disease in assoc_diseases]
                att_mappings[vertex.index].extend(assoc_disease_ids)
