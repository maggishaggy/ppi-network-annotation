import logging
from ppi_network_annotation.model.network import Network

logger = logging.getLogger(__name__)


class LabeledNetwork:
    """Mimic encapsulation of a labeled and annotated PPI network"""

    def __init__(self, network: Network):
        """Initialize the network object

        :param network: A PPI network annotated with differential gene expression and disease association.
        """
        logger.info("Initializing LabeledNetwork")

        self.graph = network.graph

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