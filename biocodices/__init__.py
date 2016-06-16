import os


if 'DISPLAY' not in os.environ:
    import matplotlib
    # Force matplotlib to not use any Xwindows backend.
    matplotlib.use('Agg')
    # This prevents matplotlib from raising an exception when running
    # on a remote server with no X. This line has to be executed before
    # importing pyplot.

__version__ = '0.3.1'
__program_name__ = 'biocodices'
software_name = '{} {}'.format(__program_name__, __version__)
