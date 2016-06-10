import matplotlib

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
# This prevents matplotlib raising an exception when running biocodices on a
# remote server with no X. This line has to be executed before importing
# pyplot, so we have to run it before any biocodices code is imported.

__version__ = '0.3'
__program_name__ = 'biocodices'
software_name = '{} {}'.format(__program_name__, __version__)

from .components import Sample, Cohort, Project
from .analyzers import AssociationTester
from .helpers.db import DB
#  from .components.dataset import Dataset
#  from .components.sample_group import SampleGroup
#  from .components.panel import Panel
#  from .helpers.config import Config
#  from .helpers.resource import Resource
