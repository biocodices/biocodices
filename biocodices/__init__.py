import matplotlib

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
# This prevents matplotlib raising an exception when running biocodices on a
# remote server with no X. This line has to be executed before importing
# pyplot, so we have to run it before any biocodices code is imported.

__version__ = '0.3'
__program_name__ = 'biocodices'
software_name = '{} {}'.format(__program_name__, __version__)
