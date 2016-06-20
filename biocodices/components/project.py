from os import mkdir
from os.path import join, expanduser, abspath, basename, isdir
from glob import glob


class Project:
    def __init__(self, base_dir):
        self.dir = abspath(expanduser(base_dir))
        self.id = basename(self.dir)
        self.data_dir = join(self.dir, 'data')
        self.results_dir = join(self.dir, 'results')

        if not isdir(self.results_dir):
            mkdir(self.results_dir)

    def __repr__(self):
        return '<{} "{}">'.format(self.__class__.__name__, self.id)

    def data_files(self, pattern=None):
        return [fn for fn in glob(join(self.data_dir, (pattern or '*')))]

    def results_files(self, pattern=None):
        return [fn for fn in glob(join(self.results_dir, (pattern or '*')))]

    def results_file(self, filename):
        return join(self.results_dir, filename)

    def data_file(self, filename):
        return join(self.data_dir, filename)
