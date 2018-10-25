"""Library for annotating a protein protein interaction network with differential gene expression."""

from ppi_network_annotation.constants import get_version  # noqa: F401
from ppi_network_annotation.model import AttributeNetwork  # noqa: F401
from ppi_network_annotation.model import FilteredNetwork  # noqa: F401
from ppi_network_annotation.model import Gene  # noqa: F401
from ppi_network_annotation.model import LabeledNetwork  # noqa: F401
from ppi_network_annotation.model import Network  # noqa: F401
from ppi_network_annotation.pipeline import generate_ppi_network, parse_dge  # noqa: F401

__version__ = '0.0.3'

__title__ = 'ppi_network_annotation'
__description__ = 'Library for annotating a protein protein interaction network with differential ' \
                  'gene expression.'
__url__ = 'https://github.com/GuiltyTargets/ppi-network-annotation'

__author__ = 'Özlem Muslu'
__email__ = 'ozlemmuslu@gmail.com'

__license__ = 'MIT License'
__copyright__ = 'Copyright (c) 2018 Özlem Muslu'
