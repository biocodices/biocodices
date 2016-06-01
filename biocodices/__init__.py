__version__ = '0.3'
software_name = 'biocodices {}'.format(__version__)

from .components import Sample, Cohort
from .analyzers import AssociationTester
from .models.db import DB
#  from .components.dataset import Dataset
#  from .components.sample_group import SampleGroup
#  from .components.panel import Panel
#  from .helpers.config import Config
#  from .helpers.resource import Resource
