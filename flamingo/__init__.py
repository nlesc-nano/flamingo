"""Library API"""

import logging

from .__version__ import __version__
from .features import atomic_features, featurizer

logging.getLogger(__name__).addHandler(logging.NullHandler())
