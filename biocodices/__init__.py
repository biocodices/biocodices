__version__ = '1.2'
software_name = 'biocodices {}'.format(__version__)

from .components.sequencing import Sequencing
from .components.sample import Sample
#  from .components.dataset import Dataset
#  from .components.sample_group import SampleGroup
#  from .components.panel import Panel
#  from .helpers.config import Config
#  from .helpers.resource import Resource
